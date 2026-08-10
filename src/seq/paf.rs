use std::{
    io::BufRead,
    path::Path,
};
use smallvec::SmallVec;
use crate::{
    seq::{
        ContigId, ContigNames,
        cigar::Cigar,
        aln::Strand,
    },
    err::{error, add_path},
    ext,
};

/// PAF file reader.
pub struct PafFile<R> {
    file: R,
    line: String,
}

impl<R: BufRead> PafFile<R> {
    pub fn new(file: R) -> Self {
        Self {
            file,
            line: String::new(),
        }
    }

    /// Use this function instead of following `std::Iterator` trait to avoid data reallocation.
    pub fn next<'a>(
        &'a mut self,
        contigs: &ContigNames,
        save_tags: bool,
    ) -> crate::Result<Option<PafParseResult<'a>>> {
        loop {
            self.line.clear();
            return match self.file.read_line(&mut self.line) {
                Ok(0) => Ok(None),
                Ok(_n) => {
                    if self.line.is_empty() || self.line.starts_with('#') {
                        continue
                    }
                    if self.line.ends_with('\n') {
                        self.line.pop();
                        if self.line.ends_with('\r') {
                            self.line.pop();
                        }
                    }
                    Some(PafEntry::parse(&self.line, contigs, save_tags)).transpose()
                }
                Err(e) => Err(crate::Error::Io(e, Vec::new())),
            }
        }
    }

    #[inline(always)]
    pub fn next_wo_tags<'a>(
        &'a mut self,
        contigs: &ContigNames,
    ) -> crate::Result<Option<PafParseResult<'a>>> {
        self.next(contigs, false)
    }

    // #[inline(always)]
    // pub fn next_with_tags<'a>(
    //     &'a mut self,
    //     contigs: &ContigNames,
    // ) -> crate::Result<Option<PafParseResult<'a>>> {
    //     self.next(contigs, true)
    // }
}

#[derive(Clone, Debug)]
pub struct PafEntry {
    query_id: ContigId,
    query_len: u32,
    query_start: u32,
    query_end: u32,
    strand: Strand,

    target_id: ContigId,
    target_len: u32,
    target_start: u32,
    target_end: u32,

    n_matches: u32,
    aln_len: u32,
    cigar: Option<Cigar>,
    // oth_tags: String,
}

#[derive(Debug)]
pub enum PafParseResult<'a> {
    Entry(PafEntry),
    UnknownContig(&'a str),
}

impl PafEntry {
    /// Tries to parse a PAF line, returns Err if parsing failed, `Ok(UnknownContig)` if contig names
    /// could not be found, or `Ok(Entry)` otherwise.
    /// If `save_tags` is false, does not save tags other than cigar.
    pub fn parse<'a>(
        line: &'a str,
        contigs: &ContigNames,
        _save_tags: bool,
    ) -> crate::Result<PafParseResult<'a>> {
        let parse_error = || error!(ParsingError, "Could not parse PAF line `{}`", line.trim());

        use PafParseResult::*;
        let split: SmallVec<[&'a str; 32]> = line.split('\t').collect();
        if split.len() < 12 {
            return Err(error!(InvalidInput, "PAF line ({}) has too few columns", line));
        }
        let query_name = split[0];
        let Some(query_id) = contigs.try_get_id(query_name) else { return Ok(UnknownContig(query_name)) };
        let target_name = split[5];
        let Some(target_id) = contigs.try_get_id(target_name) else { return Ok(UnknownContig(target_name)) };

        let query_len = split[1].parse().map_err(|_| parse_error())?;
        let query_start = split[2].parse().map_err(|_| parse_error())?;
        let query_end = split[3].parse().map_err(|_| parse_error())?;
        let strand = Strand::from_str(split[4])?;
        let target_len = split[6].parse().map_err(|_| parse_error())?;
        let target_start = split[7].parse().map_err(|_| parse_error())?;
        let target_end = split[8].parse().map_err(|_| parse_error())?;
        let n_matches = split[9].parse().map_err(|_| parse_error())?;
        let aln_len = split[10].parse().map_err(|_| parse_error())?;

        let mut entry = Self { query_id, query_len, query_start, query_end, strand,
            target_id, target_len, target_start, target_end, n_matches, aln_len,
            cigar: None,
            // oth_tags: String::new(),
        };
        for tag in &split[12..] {
            if tag.starts_with("cg:Z:") {
                entry.cigar = Some(Cigar::from_str(&tag.as_bytes()[5..])?);
            } // else if save_tags {
            //     entry.push_tag(tag);
            // }
        }
        Ok(Entry(entry))
    }

    #[inline(always)]
    pub fn query_id(&self) -> ContigId {
        self.query_id
    }

    #[inline(always)]
    pub fn query_len(&self) -> u32 {
        self.query_len
    }

    #[inline(always)]
    pub fn query_start(&self) -> u32 {
        self.query_start
    }

    #[inline(always)]
    pub fn query_end(&self) -> u32 {
        self.query_end
    }

    #[inline(always)]
    pub fn strand(&self) -> Strand {
        self.strand
    }

    #[inline(always)]
    pub fn target_id(&self) -> ContigId {
        self.target_id
    }

    #[inline(always)]
    pub fn target_len(&self) -> u32 {
        self.target_len
    }

    #[inline(always)]
    pub fn target_start(&self) -> u32 {
        self.target_start
    }

    #[inline(always)]
    pub fn target_end(&self) -> u32 {
        self.target_end
    }

    #[inline(always)]
    pub fn n_matches(&self) -> u32 {
        self.n_matches
    }

    #[inline(always)]
    pub fn aln_len(&self) -> u32 {
        self.aln_len
    }

    #[inline(always)]
    pub fn edit_distance(&self) -> Option<u32> {
        if self.aln_len == 0 { None } else { Some(self.aln_len - self.n_matches) }
    }

    #[inline(always)]
    pub fn divergence(&self) -> Option<f64> {
        self.edit_distance().map(|d| f64::from(d) / f64::from(self.aln_len))
    }

    /// Checks if the alignment fully covers both sequences and is on the forward strand.
    pub fn full_positive_alignment(&self) -> bool {
        self.strand().is_forward()
            && self.query_start() == 0 && self.query_end() == self.query_len()
            && self.target_start() == 0 && self.target_end() == self.target_len()
    }

    #[inline(always)]
    pub fn cigar(&self) -> Option<&Cigar> {
        self.cigar.as_ref()
    }

    pub fn take_cigar(&mut self) -> Option<Cigar> {
        self.cigar.take()
    }

    // pub fn other_tags(&self) -> &str {
    //     &self.oth_tags
    // }

    // pub fn push_tag(&mut self, tag: &str) {
    //     self.oth_tags.push('\t');
    //     self.oth_tags.push_str(tag);
    // }
}

/// Prune PAF file to only keep given contig IDs.
pub fn prune_paf(
    in_filename: &Path,
    out_filename: &Path,
    contigs: &ContigNames,
    keep_ids: &[ContigId],
) -> crate::Result<()> {
    let mut keep = vec![false; contigs.len()];
    for &id in keep_ids {
        keep[id.ix()] = true;
    }
    let in_paf = ext::sys::open(in_filename)?;
    let mut out_paf = ext::sys::create(out_filename)?;

    for line in in_paf.lines() {
        let line = line.map_err(add_path!(in_filename))?;
        if line.starts_with('#') {
            writeln!(out_paf, "{}", line).map_err(add_path!(out_filename))?;
            continue;
        }
        let split: SmallVec<[&str; 7]> = line.splitn(7, '\t').collect();
        if split.len() < 7 {
            return Err(error!(ParsingError, "Could not parse PAF line `{}`", line.trim()));
        }
        let keep1 = contigs.try_get_id(&split[0]).map(|i| keep[i.ix()]).unwrap_or(false);
        let keep2 = contigs.try_get_id(&split[5]).map(|i| keep[i.ix()]).unwrap_or(false);
        if keep1 && keep2 {
            writeln!(out_paf, "{}", line).map_err(add_path!(out_filename))?;
        }
    }
    Ok(())
}
