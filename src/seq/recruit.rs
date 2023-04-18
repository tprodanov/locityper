//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    collections::hash_map::Entry,
    path::Path,
    fs::File,
    io::{self, BufWriter},
    sync::mpsc::{self, Sender, Receiver},
};
use nohash::IntMap;
use flate2::write::GzEncoder;
use smallvec::{smallvec, SmallVec};
use crate::{
    Error,
    ext::{sys as sys_ext},
    seq::{
        ContigSet,
        kmers::{KmerCount, canonical_kmers, N_KMER},
        fastx::RecordExt,
    },
};

const CAPACITY: usize = 2;

/// Each k-mer is associated with a vector of loci indices.
type KmersToLoci = IntMap<u64, SmallVec<[u16; CAPACITY]>>;

/// Recruitement targets: sequences and output FASTQ files.
pub struct Targets {
    k: u8,
    out_files: Vec<BufWriter<GzEncoder<File>>>,
    kmers: KmersToLoci,
}

impl Targets {
    pub fn new<'a>(
        sets: impl Iterator<Item = (&'a ContigSet, &'a Path)>,
        max_kmer_freq: KmerCount,
    ) -> Result<Self, Error>
    {
        log::info!("Generating recruitment targets");
        let mut out_files = Vec::new();
        let mut kmers = KmersToLoci::default();
        let mut curr_kmers = Vec::new();
        let mut k = 0_u8;

        let mut discarded = 0;
        for (set_ix, (contig_set, out_path)) in sets.enumerate() {
            let kmer_counts = contig_set.kmer_counts();
            if set_ix == 0 {
                k = kmer_counts.k().try_into().unwrap();
            } else {
                assert_eq!(kmer_counts.k() as u8, k, "Different loci have different k-mer sizes!");
            }

            let set_ix_u16 = u16::try_from(set_ix).expect("Too many contig sets!");
            out_files.push(sys_ext::create_gzip(out_path)?);
            for contig_id in contig_set.contigs().ids() {
                curr_kmers.clear();
                canonical_kmers(contig_set.get_seq(contig_id), k, &mut curr_kmers);
                for (&kmer, &count) in curr_kmers.iter().zip(kmer_counts.get(contig_id)) {
                    if count <= max_kmer_freq && kmer != N_KMER {
                        match kmers.entry(kmer) {
                            Entry::Occupied(entry) => {
                                let mut loci_ixs = entry.into_mut();
                                if *loci_ixs.last().unwrap() != set_ix_u16 {
                                    loci_ixs.push(set_ix_u16);
                                }
                            }
                            Entry::Vacant(entry) => {
                                entry.insert(smallvec![set_ix_u16]);
                            }
                        }
                    } else {
                        discarded += 1;
                    }
                }
            }
        }
        log::info!("Collected {} k-mers, discarded {} k-mers", kmers.len(), discarded);
        Ok(Self { k, kmers, out_files })
    }
}

// struct Worker {
//     worker_id: u16,
//     out_files: Vec<BufWriter<File>>,
//     kmers: IntMap<u64, SmallVec<[u16; CAPACITY]>>,
//     receiver
// }

/// Record a specific single- or paired-read to one or more loci.
/// `rec_kmers` and `occurances` are buffers, filled anew for each record in order to reduce the number of allocs.
///
/// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
/// for which the number of occurances will be at least `threshold`.
fn recruit<W: io::Write>(
    record: &impl RecordExt,
    out_files: &mut [W],
    k: u8,
    kmers_to_loci: &KmersToLoci,
    rec_kmers: &mut Vec<u64>,
    occurances: &mut IntMap<u16, u16>,
    threshold: u16,
) -> io::Result<()>
{
    record.canonical_kmers(k, rec_kmers);
    occurances.clear();
    for kmer in rec_kmers.iter() {
        if let Some(loci_ixs) = kmers_to_loci.get(kmer) {
            for &locus_ix in loci_ixs.iter() {
                occurances.entry(locus_ix).and_modify(|counter| *counter += 1).or_insert(1);
            }
        }
    }
    for (&locus_ix, &count) in occurances.iter() {
        if count >= threshold {
            record.write_to(&mut out_files[usize::from(locus_ix)])?;
        }
    }
    Ok(())
}

// fn recruit_single_thread()
