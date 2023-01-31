use std::{
    cmp::{min, max},
    rc::Rc,
    fmt,
    mem::{self, MaybeUninit},
};
use htslib::bam::Record;
use statrs::distribution::Discrete;
use crate::{
    seq::{
        contigs::ContigNames,
        interv::Interval,
        cigar::Operation,
        aln::{ReadEnd, Alignment},
        compl::linguistic_complexity_k3,
    },
    math::{
        nbinom::{NBinom, CachedDistr},
        mnom::Multinomial,
    }
};

const N_COMPLEXITIES: usize = 2;

#[derive(Clone, Default)]
struct ReadSubProfile {
    matches: u16,
    mismatches: u16,
    insertions: u16,
    clipping: u16,
    deletions: u16,
}

impl fmt::Debug for ReadSubProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let qlen = self.matches + self.mismatches + self.insertions + self.clipping;
        write!(f, "qlen: {:3}, M: {:3}, X: {:3}, I: {:3}, S: {:3}, D: {:3}",
        qlen, self.matches, self.mismatches, self.insertions, self.clipping, self.deletions)
    }
}

#[derive(Clone, Default)]
struct ReadProfile([ReadSubProfile; N_COMPLEXITIES]);

impl ReadProfile {
    fn calculate(aln: &Alignment, read_compl: &[u8]) -> Self {
        let mut profiles: [ReadSubProfile; N_COMPLEXITIES] = Default::default();
        let cigar = aln.cigar();
        let (left_clip, right_clip) = cigar.true_clipping();
        let read_len = cigar.query_len();
        assert_eq!(read_len, read_compl.len() as u32, "Unexpected read length!");

        profiles[read_compl[0] as usize].clipping += u16::try_from(left_clip).unwrap();
        profiles[read_compl[read_compl.len() - 1] as usize].clipping += u16::try_from(right_clip).unwrap();

        // read_start and read_end: read boundaries after we remove clipping.
        // Need to do this, because we do not simply remove Soft clipping,
        // but also remove all operations before the first (last) match.
        let read_start = left_clip;
        let read_end = read_len - right_clip;
        let mut read_pos = 0;
        for item in cigar.iter() {
            let curr_len = item.len();
            match item.operation() {
                Operation::Match => panic!("Expected an extended CIGAR, but found M operation."),
                Operation::Equal => {
                    for i in max(read_start, read_pos)..min(read_end, read_pos + curr_len) {
                        profiles[read_compl[i as usize] as usize].matches += 1;
                    }
                    read_pos += curr_len;
                },
                Operation::Diff => {
                    // Same as in Operation::Equal.
                    for i in max(read_start, read_pos)..min(read_end, read_pos + curr_len) {
                        profiles[read_compl[i as usize] as usize].mismatches += 1;
                    }
                    read_pos += curr_len;
                },
                Operation::Ins => {
                    // No need to check read_pos + curr_len <= read_end, as insert will either be completely inside
                    // [read_start, read_end), or will be completely out.
                    if read_start <= read_pos && read_pos < read_end {
                        for i in read_pos..read_pos + curr_len {
                            profiles[read_compl[i as usize] as usize].insertions += 1;
                        }
                    }
                    read_pos += curr_len;
                },
                Operation::Soft => {
                    // Soft clipping is already accounted for.
                    read_pos += curr_len;
                },
                Operation::Del => {
                    // No need to check read_pos + curr_len <= read_end, as the deletion
                    // will either be completely inside [read_start, read_end), or will be completely out.
                    if read_start <= read_pos && read_pos < read_end {
                        profiles[read_compl[read_pos as usize] as usize].deletions += u16::try_from(curr_len).unwrap();
                    }
                },
            }
        }
        Self(profiles)
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
    insertions: u64,
    clipping: u64,
    deletions: u64,
}

impl SubProfileBuilder {
    fn update(&mut self, read_subprofile: &ReadSubProfile) {
        self.matches += u64::from(read_subprofile.matches);
        self.mismatches += u64::from(read_subprofile.mismatches);
        self.insertions += u64::from(read_subprofile.insertions);
        self.clipping += u64::from(read_subprofile.clipping);
        self.deletions += u64::from(read_subprofile.deletions);
    }

    fn sum_read_len(&self) -> u64 {
        self.matches + self.mismatches + self.insertions + self.clipping
    }

    fn get_profile(&self) -> SubErrorProfile {
        let sum_read_lenf = self.sum_read_len() as f64;
        SubErrorProfile {
            read_operation_distr: Multinomial::new(&[
                self.matches as f64 / sum_read_lenf,
                self.mismatches as f64 / sum_read_lenf,
                self.insertions as f64 / sum_read_lenf,
                self.clipping as f64 / sum_read_lenf,
                ]),
            no_del_prob: sum_read_lenf / (sum_read_lenf + self.deletions as f64),
        }
    }
}

impl fmt::Debug for SubProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sum_read_len = self.sum_read_len();
        let sum_read_lenf = sum_read_len as f64;
        write!(f, concat!("sum read length: {:8},  M: {:8} ({:.6}),  X: {:8} ({:.6}),  I: {:8} ({:.6}),  ",
                            "S: {:8} ({:.6}),  D: {:8} ({:.6})"),
            sum_read_len, self.matches, self.matches as f64 / sum_read_lenf,
            self.mismatches, self.mismatches as f64 / sum_read_lenf,
            self.insertions, self.insertions as f64 / sum_read_lenf,
            self.clipping, self.clipping as f64 / sum_read_lenf,
            self.deletions, self.deletions as f64 / (sum_read_lenf + self.deletions as f64))
    }
}

use crate::seq::aln::READ_ENDS;
type SubErrorProfiles = [[SubErrorProfile; N_COMPLEXITIES]; READ_ENDS];

/// Error profile builder.
#[derive(Default)]
struct ProfileBuilder([[SubProfileBuilder; N_COMPLEXITIES]; READ_ENDS]);

impl ProfileBuilder {
    fn update(&mut self, read_prof: &ReadProfile, read_end: ReadEnd) {
        let subprofiles = &mut self.0[read_end.ix()];
        for i in 0..N_COMPLEXITIES {
            subprofiles[i].update(&read_prof.0[i]);
        }
    }

    fn get_profiles(&self) -> SubErrorProfiles {
        // Start with uninitialized array and fill it one by one.
        // Perhaps, later there will be a simpler solution in Rust.
        let mut uninit_subprofiles: [[MaybeUninit<SubErrorProfile>; N_COMPLEXITIES]; READ_ENDS] = unsafe {
            MaybeUninit::uninit().assume_init()
        };
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                uninit_subprofiles[read_end][compl] = MaybeUninit::new(self.0[read_end][compl].get_profile());
            }
        }
        unsafe { mem::transmute::<_, SubErrorProfiles>(uninit_subprofiles) }
    }
}

impl fmt::Debug for ProfileBuilder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                writeln!(f, "Mate {},  compl {}", read_end + 1, compl)?;
                writeln!(f, "    {:?}", self.0[read_end][compl])?;
            }
        }
        Ok(())
    }
}

struct SubErrorProfile {
    /// Distribution of matches, mismatches, insertions and clippings.
    read_operation_distr: Multinomial,
    /// Probability of success (no deletion) per base-pair in the read sequence.
    /// Used in Neg.Binomial distribution.
    no_del_prob: f64,
}

impl fmt::Debug for SubErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SubErrorProfile{{ MXIS: {}, D: NB({:.6}) }}", self.read_operation_distr, self.no_del_prob)
    }
}

/// Simple structure that calculates sequence complexity in moving windows.
#[derive(Debug, Clone)]
pub struct ComplexityCalculator {
    /// Read sequence complexity is calculated on windows of corresponding size.
    window: usize,
    /// Sequence is considered *low-complexity* if complexity is not greater `thresh` (or contains Ns).
    thresh: f32,
}

impl ComplexityCalculator {
    /// Creates a new complexity calculator.
    pub fn new(window: usize, thresh: f32) -> Self {
        Self { window, thresh }
    }

    /// Calculates sequence complexity. Output size = read length, output values are in `0..N_COMPLEXITIES`.
    pub fn calculate(&self, seq: &[u8]) -> Vec<u8> {
        linguistic_complexity_k3(seq, self.window).map(|v| u8::from(!v.is_nan() && v > self.thresh)).collect()
    }
}

pub struct ErrorProfile {
    compl_calc: ComplexityCalculator,
    subprofiles: SubErrorProfiles,
}

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I>(compl_calc: ComplexityCalculator, contigs: &Rc<ContigNames>, records: I) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
    {
        log::info!("    Estimating read error profiles");
        let mut prof_builder = ProfileBuilder::default();
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            assert_ne!(record.seq_len(), 0, "Read {} does not have a sequence", String::from_utf8_lossy(record.qname()));
            let aln = Alignment::from_record(record, Rc::clone(contigs));
            let read_compl = compl_calc.calculate(&record.seq().as_bytes());
            let read_prof = ReadProfile::calculate(&aln, &read_compl);
            prof_builder.update(&read_prof, ReadEnd::from_record(record));
        }

        Self {
            compl_calc,
            subprofiles: prof_builder.get_profiles(),
        }
    }
}

impl fmt::Debug for ErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                writeln!(f, "Mate{}, compl{}:  {:?}", read_end + 1, compl, self.subprofiles[read_end][compl])?;
            }
        }
        Ok(())
    }
}