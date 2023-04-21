//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    collections::hash_map::Entry,
    path::{Path, PathBuf},
    fs::File,
    io::{self, BufWriter},
    time::Instant,
    sync::mpsc::{self, Sender, Receiver},
};
use const_format::formatcp;
use nohash::IntMap;
use smallvec::{smallvec, SmallVec};
use crate::{
    Error,
    ext::{sys as sys_ext},
    seq::{
        ContigSet,
        kmers,
        fastx::{RecordExt, FastxRead},
    },
};

const CAPACITY: usize = 4;
type MinimToLoci = IntMap<u32, SmallVec<[u16; CAPACITY]>>;

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    minimizer_k: u8,
    minimizer_n: u8,
    /// Minimizers appearing across the targets.
    minim_to_loci: MinimToLoci,
    /// Minimal number of k-mer matches per read/read-pair, needed to recruit the read to the corresponding locus.
    min_matches: u8,
    n_loci: u16,
}

impl Targets {
    pub fn new<'a>(
        sets: impl Iterator<Item = &'a ContigSet>,
        minimizer_k: u8,
        minimizer_n: u8,
        min_matches: u8,
    ) -> io::Result<Self> {
        log::info!("Generating recruitment targets");
        let mut minim_to_loci = MinimToLoci::default();
        let mut buffer = Vec::new();

        let mut n_seqs = 0;
        let mut locus_ix: u16 = 0;
        for contig_set in sets {
            for contig_id in contig_set.contigs().ids() {
                buffer.clear();
                kmers::minimizers(contig_set.get_seq(contig_id), minimizer_k, minimizer_n, &mut buffer);
                for &minimizer in buffer.iter() {
                    match minim_to_loci.entry(minimizer) {
                        Entry::Occupied(entry) => {
                            let mut loci_ixs = entry.into_mut();
                            if *loci_ixs.last().unwrap() != locus_ix {
                                loci_ixs.push(locus_ix);
                            }
                        }
                        Entry::Vacant(entry) => { entry.insert(smallvec![locus_ix]); }
                    }
                }
            }
            locus_ix = locus_ix.checked_add(1)
                .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
            n_seqs += contig_set.len();
        }
        let n_loci = locus_ix;
        log::info!("Collected {} minimizers across {} loci and {} sequences", minim_to_loci.len(), n_loci, n_seqs);
        Ok(Self { minimizer_k, minimizer_n, minim_to_loci, min_matches, n_loci })
    }

    /// Opens output fastq file for each locus.
    fn create_out_files(
        &self,
        out_filenames: impl Iterator<Item = impl AsRef<Path>>
    ) -> io::Result<Vec<impl io::Write>>
    {
        let out_files: Vec<_> = out_filenames
            .map(|path| std::fs::File::create(path.as_ref()).map(io::BufWriter::new))
            .collect::<Result<_, io::Error>>()?;
        assert_eq!(out_files.len(), usize::from(self.n_loci),
            "The number of output files does not match the number of loci");
        Ok(out_files)
    }

    /// Record a specific single- or paired-read to one or more loci.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    /// Returns true if the read was recruited anywhere.
    fn recruit_record<T: RecordExt>(
        &self,
        record: &T,
        loci_ixs: &mut Vec<u16>,
        buffer_minimizers: &mut Vec<u32>,
        buffer_matches: &mut IntMap<u16, u8>,
    )
    {
        // log::debug!("    Recruit record");
        buffer_minimizers.clear();
        record.minimizers(self.minimizer_k, self.minimizer_n, buffer_minimizers);
        for minimizer in buffer_minimizers.iter() {
            if let Some(loci_ixs) = self.minim_to_loci.get(minimizer) {
                for &locus_ix in loci_ixs.iter() {
                    buffer_matches.entry(locus_ix)
                        .and_modify(|counter| *counter = counter.saturating_add(1))
                        .or_insert(1);
                }
            }
        }
        if !buffer_matches.is_empty() {
            // log::debug!("    Matches: {:?}", buffer_matches);
            loci_ixs.clear();
            for (&locus_ix, &count) in buffer_matches.iter() {
                if count >= self.min_matches {
                    loci_ixs.push(locus_ix);
                }
            }
            // log::debug!("    Output: {:?}", loci_ixs);
            buffer_matches.clear();
        }
    }

    fn recruit_single_thread<T: RecordExt>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        out_filenames: impl Iterator<Item = impl AsRef<Path>>,
    ) -> io::Result<()>
    {
        let mut writers = self.create_out_files(out_filenames)?;
        let mut record = T::default();
        let mut read_loci = Vec::new();
        let mut buffer1 = Vec::new();
        let mut buffer2 = IntMap::default();

        let timer = Instant::now();
        let mut processed = 0_u64;
        let mut recruited = 0_u64;
        while reader.read_next(&mut record)? {
            self.recruit_record(&record, &mut read_loci, &mut buffer1, &mut buffer2);
            if !read_loci.is_empty() {
                for &locus_ix in read_loci.iter() {
                    record.write_to(&mut writers[usize::from(locus_ix)])?;
                }
                recruited += 1;
            }
            processed += 1;
            if processed % 100_000 == 0 {
                let per_read = timer.elapsed().as_secs_f64() / processed as f64;
                log::debug!("    Recruited {:11} /{:11} reads,  {:.2} us/read", recruited, processed,
                    per_read * 1e6);
            }
        }
        log::info!("Finished recruitment in {}", crate::ext::fmt::Duration(timer.elapsed()));
        Ok(())
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    pub fn recruit<T: RecordExt>(
        &self,
        reader: impl FastxRead<Record = T>,
        out_filenames: impl Iterator<Item = impl AsRef<Path>>,
        threads: u16,
    ) -> io::Result<()>
    {
        if threads <= 1 {
            self.recruit_single_thread(reader, out_filenames)
        } else {
            unimplemented!("Multi-thread recruitment is not implemented yet.")
        }
    }
}
