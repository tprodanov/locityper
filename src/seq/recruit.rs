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
        kmers::{KmerCount, canonical_kmers, N_KMER},
        fastx::{RecordExt, FastxRead},
    },
};

// fn murmur3_hash(mut key: u64) -> u64 {
//     key ^= key >> 33;
//     key *= 0xff51afd7ed558ccd;
//     key ^= key >> 33;
//     key *= 0xc4ceb9fe1a85ec53;
//     key ^= key >> 33
//     key
// }

// fn get_minimizers(&kmers: &[u64], ) ->

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    /// k-mer size.
    k: u8,
    /// k-mers appearing across the targets.
    kmers_set: IntSet<u64>,
    /// Number of contig sets.
    n_loci: u16,
    /// Minimal number of k-mer matches per read/read-pair, needed to recruit the read to the corresponding locus.
    min_matches: u16,
}

impl Targets {
    pub fn new<'a>(
        sets: impl Iterator<Item = &'a ContigSet>,
        max_kmer_freq: KmerCount,
        min_matches: u16
    ) -> io::Result<Self> {
        log::info!("Generating recruitment targets");
        let mut kmers_set = IntSet::default();
        let mut discarded_kmers = IntSet::default();
        let mut curr_kmers = Vec::new();
        let mut k = 0_u8;

        let mut set_ix: u16 = 0;
        for contig_set in sets {
            let kmer_counts = contig_set.kmer_counts();
            if set_ix == 0 {
                k = kmer_counts.k().try_into().unwrap();
            } else {
                assert_eq!(kmer_counts.k() as u8, k, "Different loci have different k-mer sizes!");
            }

            for contig_id in contig_set.contigs().ids() {
                curr_kmers.clear();
                canonical_kmers(contig_set.get_seq(contig_id), k, &mut curr_kmers);
                for (&kmer, &count) in curr_kmers.iter().zip(kmer_counts.get(contig_id)) {
                    if count <= max_kmer_freq && kmer != N_KMER {
                        kmers_set.insert(kmer);
                    } else {
                        discarded_kmers.insert(kmer);
                    }
                }
            }
            set_ix = set_ix.checked_add(1)
                .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
        }
        let n_loci = set_ix;
        log::info!("Collected {} k-mers across {} loci, discarded {} k-mers",
            kmers_set.len(), n_loci, discarded_kmers.len());
        Ok(Self { k, kmers_set, n_loci, min_matches })
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
    /// `rec_kmers` and `loci_matches` are buffers, filled anew for each record in order to reduce the number of allocs.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    /// Returns true if the read was recruited anywhere.
    fn recruit_record<T: RecordExt>(
        &self,
        record: &T,
        rec_kmers: &mut Vec<u64>,
    ) -> bool
    {
        record.canonical_kmers(self.k, rec_kmers);
        let mut matches = 0;
        for kmer in rec_kmers.iter() {
            if self.kmers_set.contains(kmer) {
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
        let mut rec_kmers = Vec::new();

        let timer = Instant::now();
        let mut processed = 0_u64;
        let mut recruited = 0_u64;
        while reader.read_next(&mut record)? {
            if self.recruit_record(&record, &mut rec_kmers) {
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
