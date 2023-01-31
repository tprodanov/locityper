use std::{
    cmp::{min, max},
    rc::Rc,
    fmt,
};
use htslib::bam::Record;
use statrs::distribution::Discrete;
use crate::{
    seq::{
        contigs::ContigNames,
        interv::Interval,
        cigar::{Operation, Cigar},
        aln::{ReadEnd, Strand, Alignment},
        compl::linguistic_complexity_k3,
    },
    math::{
        nbinom::{NBinom, CachedDistr},
        mnom::Multinomial,
    },
    bg::ser::{JsonSer, LoadError},
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

impl ReadSubProfile {
    fn read_len(&self) -> u16 {
        self.matches + self.mismatches + self.insertions + self.clipping
    }
}

impl fmt::Debug for ReadSubProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "qlen: {:3}, M: {:3}, X: {:3}, I: {:3}, S: {:3}, D: {:3}",
            self.read_len(), self.matches, self.mismatches, self.insertions, self.clipping, self.deletions)
    }
}

#[derive(Clone, Default)]
struct ReadProfile([ReadSubProfile; N_COMPLEXITIES]);

impl ReadProfile {
    fn calculate(cigar: &Cigar, read_complx: &[u8]) -> Self {
        let mut profiles: [ReadSubProfile; N_COMPLEXITIES] = Default::default();
        let (left_clip, right_clip) = cigar.true_clipping();
        let read_len = cigar.query_len();
        assert_eq!(read_len, read_complx.len() as u32, "Unexpected read length!");

        profiles[read_complx[0] as usize].clipping += u16::try_from(left_clip).unwrap();
        profiles[read_complx[read_complx.len() - 1] as usize].clipping += u16::try_from(right_clip).unwrap();

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
                        profiles[read_complx[i as usize] as usize].matches += 1;
                    }
                    read_pos += curr_len;
                },
                Operation::Diff => {
                    // Same as in Operation::Equal.
                    for i in max(read_start, read_pos)..min(read_end, read_pos + curr_len) {
                        profiles[read_complx[i as usize] as usize].mismatches += 1;
                    }
                    read_pos += curr_len;
                },
                Operation::Ins => {
                    // No need to check read_pos + curr_len <= read_end, as insert will either be completely inside
                    // [read_start, read_end), or will be completely out.
                    if read_start <= read_pos && read_pos < read_end {
                        for i in read_pos..read_pos + curr_len {
                            profiles[read_complx[i as usize] as usize].insertions += 1;
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
                        profiles[read_complx[read_pos as usize] as usize].deletions += u16::try_from(curr_len).unwrap();
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
        let mut subprofiles: SubErrorProfiles = Default::default();
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                subprofiles[read_end][compl] = self.0[read_end][compl].get_profile();
            }
        }
        subprofiles
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

#[derive(Default, Clone)]
pub struct ForwRevComplexity {
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

impl ForwRevComplexity {
    fn clear(&mut self) {
        self.forward.clear();
        self.reverse.clear();
    }

    fn clear_and_get(&mut self, strand: Strand) -> &mut Vec<u8> {
        self.forward.clear();
        self.reverse.clear();
        if strand.is_forward() {
            &mut self.forward
        } else {
            &mut self.reverse
        }
    }

    fn is_empty(&self) -> bool {
        self.forward.is_empty() && self.reverse.is_empty()
    }

    fn get(&mut self, strand: Strand) -> &[u8] {
        if strand.is_forward() {
            if self.forward.is_empty() {
                assert!(!self.reverse.is_empty(),
                    "ForwRevComplexity: both forward and reverse complexities are uninitialized!");
                self.forward.extend(self.reverse.iter().rev().cloned());
            }
            &self.forward
        } else {
            if self.reverse.is_empty() {
                assert!(!self.forward.is_empty(),
                    "ForwRevComplexity: both forward and reverse complexities are uninitialized!");
                self.reverse.extend(self.forward.iter().rev().cloned());
            }
            &self.reverse
        }
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

    /// Calculates sequence complexity, clears buffer, and write `seq.len()` elements to it.
    /// Possible output values are in `0..N_COMPLEXITIES`.
    pub fn calculate_direct(&self, seq: &[u8], buffer: &mut Vec<u8>) {
        buffer.extend(linguistic_complexity_k3(seq, self.window).map(|v| u8::from(!v.is_nan() && v > self.thresh)));
    }

    /// Lazily calculates forward/reverse sequence complexity for the record.
    /// Clears and updates `forw_rev_complx` buffer.
    pub fn calculate_forw_rev(&self, record: &Record, forw_rev_complx: &mut ForwRevComplexity) {
        let read_seq = record.seq().as_bytes();
        assert!(!read_seq.is_empty(), "Cannot calculate sequence complexity: read {} has no sequence",
            String::from_utf8_lossy(record.qname()));
        let strand = Strand::from_record(record);
        self.calculate_direct(&read_seq, forw_rev_complx.clear_and_get(strand));
    }
}

impl JsonSer for ComplexityCalculator {
    fn save(&self) -> json::JsonValue {
        json::object!{
            window: self.window,
            thresh: self.thresh,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        let window = obj["window"].as_usize().ok_or_else(|| LoadError(format!(
            "ComplexityCalculator: Failed to parse '{}': missing or incorrect 'window' field!", obj)))?;
        let thresh = obj["thresh"].as_f32().ok_or_else(|| LoadError(format!(
            "ComplexityCalculator: Failed to parse '{}': missing or incorrect 'thresh' field!", obj)))?;
        Ok(Self { window, thresh })
    }
}

#[derive(Default)]
struct SubErrorProfile {
    /// Distribution of matches, mismatches, insertions and clippings.
    read_operation_distr: Multinomial,
    /// Probability of success (no deletion) per base-pair in the read sequence.
    /// Used in Neg.Binomial distribution.
    no_del_prob: f64,
}

impl SubErrorProfile {
    fn ln_prob(&self, subprofile: &ReadSubProfile) -> f64 {
        let read_len = subprofile.read_len();
        if read_len == 0 {
            0.0
        } else {
            NBinom::new(f64::from(read_len), self.no_del_prob).ln_pmf_f64(f64::from(subprofile.deletions))
                + self.read_operation_distr.ln_pmf(
                    &[subprofile.matches, subprofile.mismatches, subprofile.insertions, subprofile.clipping])
        }
    }
}

impl fmt::Debug for SubErrorProfile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SubErrorProfile{{ MXIS: {}, D: NB({:.6}) }}", self.read_operation_distr, self.no_del_prob)
    }
}

impl JsonSer for SubErrorProfile {
    fn save(&self) -> json::JsonValue {
        json::object!{
            op_distr: self.read_operation_distr.save(),
            no_del_prob: self.no_del_prob,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if !obj.has_key("op_distr") {
            return Err(LoadError(format!(
                "SubErrorProfile: Failed to parse '{}': missing 'op_distr' field!", obj)))?;
        }
        let read_operation_distr = Multinomial::load(&obj["op_distr"])?;
        let no_del_prob = obj["no_del_prob"].as_f64().ok_or_else(|| LoadError(format!(
            "SubErrorProfile: Failed to parse '{}': missing or incorrect 'no_del_prob' field!", obj)))?;
        Ok(Self { read_operation_distr, no_del_prob })
    }
}

pub struct ErrorProfile {
    compl_calc: ComplexityCalculator,
    subprofiles: SubErrorProfiles,
}

impl ErrorProfile {
    /// Create error profile from the iterator over records.
    pub fn estimate<'a, I>(compl_calc: ComplexityCalculator, records: I) -> ErrorProfile
    where I: Iterator<Item = &'a Record>,
    {
        log::info!("    Estimating read error profiles");
        let mut prof_builder = ProfileBuilder::default();
        let mut read_complx = Vec::new();
        for record in records {
            // Unmapped or secondary.
            if (record.flags() & 3844) != 0 {
                continue;
            }
            assert_ne!(record.seq_len(), 0, "Read {} does not have a sequence", String::from_utf8_lossy(record.qname()));
            let cigar = Cigar::infer_ext_cigar(record, ());
            read_complx.clear();
            compl_calc.calculate_direct(&record.seq().as_bytes(), &mut read_complx);
            let read_prof = ReadProfile::calculate(&cigar, &read_complx);
            prof_builder.update(&read_prof, ReadEnd::from_record(record));
        }

        Self {
            compl_calc,
            subprofiles: prof_builder.get_profiles(),
        }
    }

    /// Returns complexity calculator.
    pub fn complexity_calculator(&self) -> &ComplexityCalculator {
        &self.compl_calc
    }

    /// Returns alignment ln-probability.
    pub fn ln_prob(&self, aln: &Alignment, read_end: ReadEnd, read_complx: &mut ForwRevComplexity) -> f64 {
        let read_prof = ReadProfile::calculate(aln.cigar(), read_complx.get(aln.strand()));
        self.subprofiles[read_end.ix()].iter()
            .zip(&read_prof.0)
            .map(|(err_subprofile, read_subprofile)| err_subprofile.ln_prob(read_subprofile))
            .sum()
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

impl JsonSer for ErrorProfile {
    fn save(&self) -> json::JsonValue {
        let mut subprofiles = Vec::with_capacity(READ_ENDS * N_COMPLEXITIES);
        for read_end in 0..READ_ENDS {
            for compl in 0..N_COMPLEXITIES {
                subprofiles.push(self.subprofiles[read_end][compl].save());
            }
        }
        json::object!{
            compl_calc: self.compl_calc.save(),
            subprofiles: subprofiles,
        }
    }

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError> {
        if !obj.has_key("compl_calc") || !obj.has_key("subprofiles") {
            return Err(LoadError(format!(
                "ErrorProfile: Failed to parse '{}': missing 'compl_calc' or 'subprofiles' field!", obj)))?;
        }
        let compl_calc = ComplexityCalculator::load(&obj["compl_calc"])?;

        let mut subprofiles: SubErrorProfiles = Default::default();
        let key = "subprofiles";
        if let json::JsonValue::Array(arr) = &obj[key] {
            if arr.len() != READ_ENDS * N_COMPLEXITIES {
                return Err(LoadError(format!("Failed to parse '{}': incorrect number of elements in array '{}'",
                    obj, key)));
            }
            for (i, obj) in arr.iter().enumerate() {
                subprofiles[i / READ_ENDS][i % READ_ENDS] = SubErrorProfile::load(obj)?;
            }
        } else {
            return Err(LoadError(format!("Failed to parse '{}': missing or incorrect array '{}'", obj, key)));
        }
        Ok(Self { compl_calc, subprofiles })
    }
}