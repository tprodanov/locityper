use std::{
    cmp::{min, max},
    rc::Rc,
    fmt,
};
use htslib::bam::Record;
use intmap::{IntMap, Entry};
use crate::{
    seq::{
        contigs::ContigNames,
        interv::Interval,
        cigar::Operation,
        aln::{ReadEnd, Alignment},
        compl::linguistic_complexity_k3,
    },
};

const N_COMPLEXITIES: usize = 2;

#[derive(Clone, Default)]
struct ReadSubProfile {
    matches: u16,
    mismatches: u16,
    deletions: u16,
    insertions: Vec<u16>,
    left_clip: Option<u16>,
    right_clip: Option<u16>,
}

impl ReadSubProfile {
    fn clear(&mut self) {
        self.matches = 0;
        self.mismatches = 0;
        self.deletions = 0;
        self.insertions.clear();
        self.left_clip = None;
        self.right_clip = None;
    }

    /// Sum number of reference positions, covered by this profile (matches + mismatches + deletions).
    fn ref_positions(&self) -> u16 {
        self.matches + self.mismatches + self.deletions
    }
}

impl fmt::Debug for ReadSubProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let obs = self.ref_positions();
        write!(f, "{:3} obs.  M: {:3}, X: {:3}, D: {:3}, I: {:?},  S {:>2} {:>2}",
            obs, self.matches, self.mismatches, self.deletions, self.insertions,
            self.left_clip.map(|c| c.to_string()).unwrap_or_else(|| "?".to_string()),
            self.right_clip.map(|c| c.to_string()).unwrap_or_else(|| "?".to_string()))
    }
}

#[derive(Clone, Default)]
struct ReadProfile([ReadSubProfile; N_COMPLEXITIES]);

impl ReadProfile {
    fn clear(&mut self) {
        self.0.iter_mut().for_each(ReadSubProfile::clear);
    }

    fn fill(&mut self, aln: &Alignment, compl_interval: &Interval, complexities: &[u8]) {
        let ref_interval = aln.ref_interval();
        assert!(ref_interval.overlaps(compl_interval),
            "Cannot calculate read profile: read alignment ({}) does not overlap with the reference interval ({}).",
            ref_interval, compl_interval);
        let cigar = aln.cigar();
        let (left_clip, right_clip) = cigar.true_clipping();

        let ref_start = ref_interval.start();
        let ref_end = ref_interval.end();
        let compl_start = compl_interval.start();
        let compl_end = compl_interval.end();
        if ref_start >= compl_start {
            self.0[complexities[(ref_start - compl_start) as usize] as usize].left_clip =
                Some(u16::try_from(left_clip).unwrap());
        }
        if ref_end <= compl_end {
            self.0[complexities[(ref_end - compl_start - 1) as usize] as usize].right_clip =
                Some(u16::try_from(right_clip).unwrap());
        }

        let read_start = left_clip;
        let read_end = cigar.query_len() - right_clip;
        let mut read_pos = 0;
        let mut ref_pos = ref_start;
        for item in cigar.iter() {
            let curr_len = item.len();
            match item.operation() {
                Operation::Match => panic!("Expected an extended CIGAR, but found M operation."),
                Operation::Equal => {
                    // `i in 0..curr_len` that satisfies two conditions:
                    //   (a)  compl_start <=  ref_pos + i < compl_end,
                    //   (b)  read_start  <= read_pos + i < read_end.
                    let lower = max(compl_start.saturating_sub(ref_pos), read_start.saturating_sub(read_pos));
                    let upper = min(curr_len,
                        min(compl_end.saturating_sub(ref_pos), read_end.saturating_sub(read_pos)));
                    for i in lower..upper {
                        self.0[complexities[(ref_pos + i - compl_start) as usize] as usize].matches += 1;
                    }
                    read_pos += curr_len;
                    ref_pos += curr_len;
                },
                Operation::Diff => {
                    // Same as in Operation::Equal.
                    let lower = max(compl_start.saturating_sub(ref_pos), read_start.saturating_sub(read_pos));
                    let upper = min(curr_len,
                        min(compl_end.saturating_sub(ref_pos), read_end.saturating_sub(read_pos)));
                    for i in lower..upper {
                        self.0[complexities[(ref_pos + i - compl_start) as usize] as usize].mismatches += 1;
                    }
                    read_pos += curr_len;
                    ref_pos += curr_len;
                },
                Operation::Ins => {
                    if compl_start <= ref_pos && ref_pos < compl_end {
                        self.0[complexities[(ref_pos - compl_start) as usize] as usize]
                            .insertions.push(u16::try_from(curr_len).unwrap());
                    }
                    read_pos += curr_len;
                },
                Operation::Soft => {
                    read_pos += curr_len;
                },
                Operation::Del => {
                    if read_start <= read_pos && read_pos < read_end {
                        let lower = ref_pos.saturating_sub(compl_start);
                        let upper = min(ref_pos + curr_len, compl_end) - compl_start;
                        for i in lower..upper {
                            self.0[complexities[i as usize] as usize].deletions += 1;
                        }
                    }
                    ref_pos += curr_len;
                },
            }
        }

    }
}

impl fmt::Debug for ReadProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, sub_profile) in self.0.iter().enumerate() {
            writeln!(f, "    compl {}:  {:?}", i, sub_profile)?;
        }
        Ok(())
    }
}

#[derive(Clone, Default)]
struct SubProfileBuilder {
    matches: u64,
    mismatches: u64,
    deletions: u64,
    sum_insertions: u64,
    n_insertions: u32,
}

impl SubProfileBuilder {
    fn update(&mut self, read_subprofile: &ReadSubProfile) {
        self.matches += u64::from(read_subprofile.matches);
        self.mismatches += u64::from(read_subprofile.mismatches);
        self.deletions += u64::from(read_subprofile.deletions);
        for &ins in read_subprofile.insertions.iter() {
            self.n_insertions += 1;
            self.sum_insertions += u64::from(ins);
        }
    }
}

impl fmt::Debug for SubProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let obs = self.matches + self.mismatches + self.deletions;
        let obsf = obs as f64;
        write!(f, "{:8} obs,  M: {:8} ({:.6}),  X: {:8} ({:.6}),  D: {:8} ({:.6}),  I: {:6}/{:6} ({:.6})",
            obs, self.matches, self.matches as f64 / obsf, self.mismatches, self.mismatches as f64 / obsf,
            self.deletions, self.deletions as f64 / obsf,
            self.sum_insertions, self.n_insertions, self.sum_insertions as f64 / self.n_insertions as f64)
    }
}

use crate::seq::aln::READ_ENDS;

/// Construct error profiles.
struct ProfileBuilder {
    /// Reference sequence interval.
    interval: Interval,
    /// Sequence complexity around each nucleotide of the input interval.
    /// Stores values in range 0..N_COMPLEXITIES.
    complexities: Vec<u8>,
    read_prof: ReadProfile,
    subprofiles: [[SubProfileBuilder; N_COMPLEXITIES]; READ_ENDS],
}

impl ProfileBuilder {
    fn new(interval: &Interval, interval_seq: &[u8], window: usize, compl_thresh: f32) -> Self {
        assert_eq!(interval.len(), interval_seq.len() as u32);
        let complexities = linguistic_complexity_k3(interval_seq, window)
            .map(|c| u8::from(!c.is_nan() && c > compl_thresh))
            .collect();
        println!("Profile builder on interval {}", interval);
        Self {
            interval: interval.clone(),
            complexities,
            read_prof: ReadProfile::default(),
            subprofiles: Default::default(),
        }
    }

    fn update(&mut self, aln: &Alignment, read_end: ReadEnd) {
        self.read_prof.clear();
        ReadProfile::fill(&mut self.read_prof, aln, &self.interval, &self.complexities);
        let subprofiles = &mut self.subprofiles[read_end.ix()];
        for i in 0..N_COMPLEXITIES {
            subprofiles[i].update(&self.read_prof.0[i]);
        }
    }
}

impl fmt::Debug for ProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                writeln!(f, "Mate {},  compl {}", read_end + 1, compl)?;
                writeln!(f, "    {:?}", self.subprofiles[read_end][compl])?;
            }
        }
        Ok(())
    }
}

pub struct ErrorProfile;

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I>(
        records: I,
        interval: &Interval,
        interval_seq: &[u8],
        compl_window: usize,
        compl_thresh: f32,
    ) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
    {
        log::info!("    Estimating read error profiles");
        let mut prof_builder = ProfileBuilder::new(interval, interval_seq, compl_window, compl_thresh);
        let contigs = interval.contigs();
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            let aln = Alignment::from_record(&record, Rc::clone(contigs));
            prof_builder.update(&aln, ReadEnd::from_record(&record));
        }
        println!("{:?}", prof_builder);
        ErrorProfile
    }
}
