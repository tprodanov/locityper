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
use nohash::IntSet;
use crate::{
    Error,
    ext::{sys as sys_ext},
    seq::{
        ContigSet,
        kmers,
        fastx::{RecordExt, FastxRead},
    },
};

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    minimizer_k: u8,
    minimizer_n: u8,
    /// Minimizers appearing across the targets.
    minimizers: IntSet<u32>,
    /// Minimal number of k-mer matches per read/read-pair, needed to recruit the read to the corresponding locus.
    min_matches: u16,
}

impl Targets {
    pub fn new<'a>(
        seqs: impl Iterator<Item = &'a [u8]>,
        minimizer_k: u8,
        minimizer_n: u8,
        min_matches: u16
    ) -> io::Result<Self> {
        log::info!("Generating recruitment targets");
        let mut minimizers = IntSet::default();
        let mut buffer = Vec::new();

        let mut n_seqs = 0;
        for seq in seqs {
            n_seqs += 1;
            buffer.clear();
            kmers::minimizers(seq, minimizer_k, minimizer_n, &mut buffer);
            minimizers.extend(buffer.iter().copied());
        }
        log::info!("Collected {} minimizers across {} sequences", minimizers.len(), n_seqs);
        Ok(Self { minimizer_k, minimizer_n, minimizers, min_matches })
    }

    // /// Opens output fastq file for each locus.
    // fn create_out_files(
    //     &self,
    //     out_filenames: impl Iterator<Item = impl AsRef<Path>>
    // ) -> io::Result<Vec<impl io::Write>>
    // {
    //     let out_files: Vec<_> = out_filenames
    //         .map(|path| sys_ext::create_gzip(path.as_ref()))
    //         .collect::<Result<_, io::Error>>()?;
    //     assert_eq!(out_files.len(), usize::from(self.n_loci),
    //         "The number of output files does not match the number of loci");
    //     Ok(out_files)
    // }

    /// Record a specific single- or paired-read to one or more loci.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    /// Returns true if the read was recruited anywhere.
    fn recruit_record<T: RecordExt>(
        &self,
        record: &T,
        buffer: &mut Vec<u32>,
    ) -> bool
    {
        buffer.clear();
        record.minimizers(self.minimizer_k, self.minimizer_n, buffer);
        let mut matches = 0;
        for val in buffer.iter() {
            if self.minimizers.contains(val) {
                matches += 1;
                if matches >= self.min_matches {
                    return true;
                }
            }
        }
        false
    }

    fn recruit_single_thread<T: RecordExt>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        out_filename: impl AsRef<Path>,
    ) -> io::Result<()>
    {
        let mut writer = io::BufWriter::new(std::fs::File::create(out_filename)?);
        let mut record = T::default();
        let mut buffer = Vec::new();

        let timer = Instant::now();
        let mut processed = 0_u64;
        let mut recruited = 0_u64;
        while reader.read_next(&mut record)? {
            if self.recruit_record(&record, &mut buffer) {
                record.write_to(&mut writer)?;
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
        out_filename: impl AsRef<Path>,
        threads: u16,
    ) -> io::Result<()>
    {
        if threads <= 1 {
            self.recruit_single_thread(reader, out_filename)
        } else {
            unimplemented!("Multi-thread recruitment is not implemented yet.")
        }
    }
}
