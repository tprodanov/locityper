//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    collections::hash_map::Entry,
    path::{Path, PathBuf},
    fs::File,
    io::{self, BufWriter},
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
        kmers::{KmerCount, canonical_kmers, N_KMER},
        fastx::{RecordExt, FastxRead},
    },
};

const CAPACITY: usize = 2;

/// Each k-mer is associated with a vector of loci indices.
type KmersToLoci = IntMap<u64, SmallVec<[u16; CAPACITY]>>;

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    /// k-mer size.
    k: u8,
    /// IntMap k-mer -> indices of contig sets, where the k-mer appears.
    kmers_to_loci: KmersToLoci,
    /// Number of contig sets.
    n_sets: u16,
    /// Minimal number of k-mer matches per read/read-pair, needed to recruit the read to the corresponding locus.
    min_matches: u16,
}

impl Targets {
    pub fn new(sets: &[ContigSet], max_kmer_freq: KmerCount, min_matches: u16) -> Result<Self, Error> {
        log::info!("Generating recruitment targets");
        let mut kmers_to_loci = KmersToLoci::default();
        let mut curr_kmers = Vec::new();
        let mut k = 0_u8;
        let n_sets = u16::try_from(sets.len())
            .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));

        let mut discarded = 0;
        for (set_ix, contig_set) in sets.iter().enumerate() {
            let set_ix = set_ix as u16;
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
                        match kmers_to_loci.entry(kmer) {
                            Entry::Occupied(entry) => {
                                let mut loci_ixs = entry.into_mut();
                                if *loci_ixs.last().unwrap() != set_ix {
                                    loci_ixs.push(set_ix);
                                }
                            }
                            Entry::Vacant(entry) => {
                                entry.insert(smallvec![set_ix]);
                            }
                        }
                    } else {
                        discarded += 1;
                    }
                }
            }
        }
        log::info!("Collected {} k-mers, discarded {} k-mers", kmers_to_loci.len(), discarded);
        Ok(Self { k, kmers_to_loci, n_sets, min_matches })
    }

    /// Opens output fastq file for each locus.
    fn create_out_files(&self, out_filenames: &[PathBuf]) -> io::Result<Vec<impl io::Write>> {
        assert_eq!(out_filenames.len(), usize::from(self.n_sets),
            "The number of output files does not match the number of loci");
        out_filenames.iter()
            .map(|path| sys_ext::create_gzip(path))
            .collect::<Result<_, io::Error>>()
    }

    /// Record a specific single- or paired-read to one or more loci.
    /// `rec_kmers` and `loci_matches` are buffers, filled anew for each record in order to reduce the number of allocs.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    fn recruit_record<T: RecordExt, W: io::Write>(
        &self,
        record: &T,
        out_files: &mut [W],
        rec_kmers: &mut Vec<u64>,
        loci_matches: &mut IntMap<u16, u16>,
    ) -> io::Result<()>
    {
        record.canonical_kmers(self.k, rec_kmers);
        loci_matches.clear();
        for kmer in rec_kmers.iter() {
            if let Some(loci_ixs) = self.kmers_to_loci.get(kmer) {
                for &locus_ix in loci_ixs.iter() {
                    loci_matches.entry(locus_ix).and_modify(|counter| *counter += 1).or_insert(1);
                }
            }
        }
        for (&locus_ix, &count) in loci_matches.iter() {
            if count >= self.min_matches {
                record.write_to(&mut out_files[usize::from(locus_ix)])?;
            }
        }
        Ok(())
    }

    fn recruit_single_thread<T: RecordExt>(
        &self,
        reader: &mut impl FastxRead<Record = T>,
        out_dirs: &[PathBuf],
    ) -> io::Result<()>
    {
        let mut out_files = self.create_out_files(out_dirs)?;
        let mut record = T::default();
        let mut rec_kmers = Vec::new();
        let mut loci_matches = IntMap::default();

        while reader.read_next(&mut record)? {
            self.recruit_record(&record, &mut out_files, &mut rec_kmers, &mut loci_matches)?;
        }
        Ok(())
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    pub fn recruit<T: RecordExt>(
        &self,
        reader: &mut impl FastxRead<Record = T>,
        out_dirs: &[PathBuf],
        threads: u16,
    ) -> io::Result<()>
    {
        if threads <= 1 {
            self.recruit_single_thread(reader, out_dirs)
        } else {
            unimplemented!("Multi-thread recruitment is not implemented yet.")
        }
    }
}
