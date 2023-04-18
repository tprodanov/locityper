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
        let mut kmers_to_loci = KmersToLoci::default();
        let mut curr_kmers = Vec::new();
        let mut k = 0_u8;

        let mut discarded = 0;
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
            set_ix = set_ix.checked_add(1)
                .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
        }
        let n_loci = set_ix;
        log::info!("Collected {} k-mers across {} loci, discarded {} k-mers", kmers_to_loci.len(), n_loci, discarded);
        Ok(Self { k, kmers_to_loci, n_loci, min_matches })
    }

    /// Opens output fastq file for each locus.
    fn create_out_files(
        &self,
        out_filenames: impl Iterator<Item = impl AsRef<Path>>
    ) -> io::Result<Vec<impl io::Write>>
    {
        let out_files: Vec<_> = out_filenames
            .map(|path| sys_ext::create_gzip(path.as_ref()))
            .collect::<Result<_, io::Error>>()?;
        assert_eq!(out_files.len(), usize::from(self.n_loci),
            "The number of output files does not match the number of loci");
        Ok(out_files)
    }

    /// Record a specific single- or paired-read to one or more loci.
    /// `rec_kmers` and `loci_matches` are buffers, filled anew for each record in order to reduce the number of allocs.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    /// Returns true if the read was recruited anywhere.
    fn recruit_record<T: RecordExt, W: io::Write>(
        &self,
        record: &T,
        out_files: &mut [W],
        rec_kmers: &mut Vec<u64>,
        loci_matches: &mut IntMap<u16, u16>,
    ) -> io::Result<bool>
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
        let mut recruited = false;
        for (&locus_ix, &count) in loci_matches.iter() {
            if count >= self.min_matches {
                recruited = true;
                record.write_to(&mut out_files[usize::from(locus_ix)])?;
            }
        }
        Ok(recruited)
    }

    fn recruit_single_thread<T: RecordExt>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        out_filenames: impl Iterator<Item = impl AsRef<Path>>,
    ) -> io::Result<()>
    {
        let mut out_files = self.create_out_files(out_filenames)?;
        let mut record = T::default();
        let mut rec_kmers = Vec::new();
        let mut loci_matches = IntMap::default();

        let mut processed = 0_u64;
        let mut recruited = 0_u64;
        while reader.read_next(&mut record)? {
            recruited += u64::from(self.recruit_record(&record, &mut out_files, &mut rec_kmers, &mut loci_matches)?);
            processed += 1;
            if processed % 100_000 == 0 {
                log::debug!("    Recruited {:11} /{:11} reads", recruited, processed);
            }
        }
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
