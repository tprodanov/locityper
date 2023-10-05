//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    io, thread,
    cmp::max,
    collections::hash_map::Entry,
    time::{Instant, Duration},
    sync::mpsc::{self, Sender, Receiver, TryRecvError},
};
use const_format::formatcp;
use nohash::IntMap;
use smallvec::{smallvec, SmallVec};
use crate::{
    err::{Error, validate_param, add_path},
    seq::{
        kmers::{self, Kmer},
        fastx::{self, FastxRead},
        fastx::UPDATE_SECS,
    },
};

pub type Minimizer = u64;

#[derive(Clone)]
pub struct Params {
    /// Recruit reads using k-mers of size `minimizer_k` that have the minimial hash across `minimizer_w`
    /// consecutive k-mers.
    /// Default: 15 and 10.
    pub minimizer_k: u8,
    pub minimizer_w: u8,
    /// Recruit reads that have at least this number of matching minimizers with one of the targets. Default: 2.
    pub matches_frac: f32,
    /// Pass reads between threads in chunks of this size. Default: 10000.
    pub chunk_size: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            minimizer_k: 11,
            minimizer_w: 5,
            matches_frac: 0.5,
            chunk_size: 10000,
        }
    }
}

impl Params {
    pub fn validate(&self) -> Result<(), Error> {
        validate_param!(0 < self.minimizer_k && self.minimizer_k <= Minimizer::MAX_KMER_SIZE,
            "Minimizer kmer-size must be within [1, {}]", Minimizer::MAX_KMER_SIZE);
        validate_param!(1 < self.minimizer_w && self.minimizer_w <= kmers::MAX_MINIMIZER_W,
            "Minimizer window-size must be within [2, {}]", kmers::MAX_MINIMIZER_W);
        validate_param!(self.matches_frac > 0.0 && self.matches_frac <= 1.0,
            "Matches ratio ({:.5}) cannot be zero, and cannot be over 1", self.matches_frac);
        validate_param!(self.chunk_size != 0, "Chunk size cannot be zero");
        Ok(())
    }

    pub fn set_matches_frac(&mut self, matches_frac: f32) -> Result<(), Error> {
        validate_param!(matches_frac > 0.0 && matches_frac <= 1.0,
            "Matches ratio ({:.5}) cannot be zero, and cannot be over 1", matches_frac);
        self.matches_frac = matches_frac;
        Ok(())
    }
}

/// Recruitment statistics: how long did it take, how many reads processed and how many recruited.
struct Stats {
    timer: Instant,
    recruited: u64,
    processed: u64,
    /// Last log message was at this duration since start.
    last_msg: Duration,
}

impl Stats {
    fn new() -> Self {
        Self {
            timer: Instant::now(),
            recruited: 0,
            processed: 0,
            last_msg: Duration::default(),
        }
    }

    /// Prints log message if enough time has passed.
    fn timed_print_log(&mut self) {
        let elapsed = self.timer.elapsed();
        if (elapsed - self.last_msg).as_secs() >= UPDATE_SECS {
            self.print_log_always(elapsed);
        }
    }

    fn print_log_always(&mut self, elapsed: Duration) {
        let processed = self.processed as f64;
        let speed = 1e-3 * processed / elapsed.as_secs_f64();
        log::debug!("    Recruited {:11} /{:8.2}M reads, {:4.0}k reads/s", self.recruited, 1e-6 * processed, speed);
        self.last_msg = elapsed;
    }

    fn finish(&mut self) {
        let elapsed = self.timer.elapsed();
        self.print_log_always(elapsed);
        log::info!("Finished recruitment in {}", crate::ext::fmt::Duration(elapsed));
    }
}

/// Trait-extension over single/paired reads.
pub trait RecruitableRecord : fastx::WritableRecord + Send + 'static {
    fn recruit(&self,
        targets: &Targets,
        answer: &mut Answer,
        rec_minims: &mut Vec<Minimizer>,
        matches: &mut IntMap<u16, u32>,
    );
}

impl<T: fastx::SingleRecord + fastx::WritableRecord + Send + 'static> RecruitableRecord for T {
    #[inline]
    fn recruit(&self,
        targets: &Targets,
        answer: &mut Answer,
        rec_minims: &mut Vec<Minimizer>,
        matches: &mut IntMap<u16, u32>,
    ) {
        targets.recruit_single_end_record(self.seq(), answer, rec_minims, matches);
    }
}

impl<T: fastx::SingleRecord + fastx::WritableRecord + Send + 'static> RecruitableRecord for [T; 2] {
    #[inline]
    fn recruit(&self,
        targets: &Targets,
        answer: &mut Answer,
        rec_minims: &mut Vec<Minimizer>,
        matches: &mut IntMap<u16, u32>,
    ) {
        targets.recruit_paired_end_record(self[0].seq(), self[1].seq(), answer, rec_minims, matches);
    }
}

const CAPACITY: usize = 4;
/// Key: minimizer, value: vector of loci indices, where the minimizer appears.
type MinimToLoci = IntMap<Minimizer, SmallVec<[u16; CAPACITY]>>;
/// Vector of loci indices, to which the read was recruited.
type Answer = Vec<u16>;

/// Target builder. Can be converted to targets using `finalize()`.
pub struct TargetBuilder {
    params: Params,
    locus_ix: u16,
    total_seqs: u32,
    minim_to_loci: MinimToLoci,
    buffer: Vec<Minimizer>,
}

impl TargetBuilder {
    pub fn new(params: Params) -> Self {
        Self {
            params,
            locus_ix: 0,
            total_seqs: 0,
            minim_to_loci: Default::default(),
            buffer: Vec::new(),
        }
    }

    /// Add set of locus alleles.
    pub fn add<'a>(&mut self, seqs: impl Iterator<Item = &'a [u8]>) {
        for seq in seqs {
            self.buffer.clear();
            kmers::minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, &mut self.buffer);
            for &minimizer in self.buffer.iter() {
                match self.minim_to_loci.entry(minimizer) {
                    Entry::Occupied(entry) => {
                        let loci_ixs = entry.into_mut();
                        if *loci_ixs.last().unwrap() != self.locus_ix {
                            loci_ixs.push(self.locus_ix);
                        }
                    }
                    Entry::Vacant(entry) => { entry.insert(smallvec![self.locus_ix]); }
                }
            }
            self.total_seqs += 1;
        }
        self.locus_ix = self.locus_ix.checked_add(1)
            .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
    }

    /// Finalize targets construction.
    /// Remove top `discard_minim` fraction of minimizers.
    pub fn finalize(self) -> Targets {
        let n_loci = self.locus_ix;
        let total_minims = self.minim_to_loci.len();
        log::info!("Collected {} minimizers across {} loci and {} sequences", total_minims, n_loci, self.total_seqs);
        assert!(total_minims > 0, "No minimizers for recruitment");
        Targets {
            params: self.params,
            n_loci,
            minim_to_loci: self.minim_to_loci,
        }
    }
}

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    params: Params,
    n_loci: u16,
    /// Minimizers appearing across the targets.
    minim_to_loci: MinimToLoci,
}

impl Targets {
    /// Record one single-end read to one or more loci.
    /// The read is recruited to all loci,
    /// for which the ratio of loci_matches will be at least `matches_frac`.
    fn recruit_single_end_record(
        &self,
        seq: &[u8],
        answer: &mut Answer,
        rec_minims: &mut Vec<Minimizer>,
        matches: &mut IntMap<u16, u32>,
    ) {
        matches.clear();
        rec_minims.clear();
        kmers::minimizers(seq, self.params.minimizer_k, self.params.minimizer_w, rec_minims);
        for minimizer in rec_minims.iter() {
            if let Some(loci_ixs) = self.minim_to_loci.get(minimizer) {
                for &locus_ix in loci_ixs.iter() {
                    *matches.entry(locus_ix).or_default() += 1;
                }
            }
        }

        answer.clear();
        const FLOAT_U32_MAX: f32 = u32::MAX as f32;
        let thresh = max(1, (rec_minims.len() as f32 * self.params.matches_frac).ceil().min(FLOAT_U32_MAX) as u32);
        for (&locus_ix, &count) in matches.iter() {
            if count >= thresh {
                answer.push(locus_ix);
            }
        }
    }

    /// Record one paired-end read to one or more loci.
    /// The read is recruited to all loci,
    /// for which the ratio of loci_matches will be at least `matches_frac` for both mates.
    fn recruit_paired_end_record(
        &self,
        seq1: &[u8],
        seq2: &[u8],
        answer: &mut Answer,
        rec_minims: &mut Vec<Minimizer>,
        matches: &mut IntMap<u16, u32>,
    ) {
        matches.clear();
        // First mate.
        rec_minims.clear();
        kmers::minimizers(seq1, self.params.minimizer_k, self.params.minimizer_w, rec_minims);
        let total1 = rec_minims.len();
        for minimizer in rec_minims.iter() {
            if let Some(loci_ixs) = self.minim_to_loci.get(minimizer) {
                for &locus_ix in loci_ixs.iter() {
                    let count = matches.entry(locus_ix).or_default();
                    let counts = unsafe { std::mem::transmute::<&mut u32, &mut [u16; 2]>(count) };
                    counts[0] = counts[0].saturating_add(1);
                }
            }
        }

        // Second mate.
        rec_minims.clear();
        kmers::minimizers(seq2, self.params.minimizer_k, self.params.minimizer_w, rec_minims);
        let total2 = rec_minims.len();
        for minimizer in rec_minims.iter() {
            if let Some(loci_ixs) = self.minim_to_loci.get(minimizer) {
                for locus_ix in loci_ixs.iter() {
                    // No reason to insert new loci if they did not match the first read end.
                    if let Some(count) = matches.get_mut(locus_ix) {
                        let counts = unsafe { std::mem::transmute::<&mut u32, &mut [u16; 2]>(count) };
                        counts[1] = counts[1].saturating_add(1);
                    }
                }
            }
        }

        answer.clear();
        const FLOAT_U16_MAX: f32 = u16::MAX as f32;
        let thresh1 = max(1, (total1 as f32 * self.params.matches_frac).ceil().min(FLOAT_U16_MAX) as u16);
        let thresh2 = max(1, (total2 as f32 * self.params.matches_frac).ceil().min(FLOAT_U16_MAX) as u16);
        for (&locus_ix, &count) in matches.iter() {
            let [count1, count2] = unsafe { std::mem::transmute::<u32, [u16; 2]>(count) };
            if count1 >= thresh1 && count2 >= thresh2 {
                answer.push(locus_ix);
            }
        }
    }

    fn recruit_single_thread<T: RecruitableRecord>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        writers: &mut [impl io::Write],
    ) -> Result<(), Error>
    {
        let mut record = T::default();
        let mut answer = Answer::new();
        let mut buffer1 = Vec::new();
        let mut buffer2 = IntMap::default();

        let mut stats = Stats::new();
        while reader.read_next(&mut record)? {
            record.recruit(self, &mut answer, &mut buffer1, &mut buffer2);
            for &locus_ix in answer.iter() {
                record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
            }
            stats.recruited += u64::from(!answer.is_empty());
            stats.processed += 1;
            if stats.processed % 10000 == 0 {
                stats.timed_print_log();
            }
        }
        stats.finish();
        Ok(())
    }

    fn recruit_multi_thread<T: RecruitableRecord>(
        &self,
        reader: impl FastxRead<Record = T>,
        writers: Vec<impl io::Write>,
        threads: u16,
    ) -> Result<(), Error> {
        let n_workers = usize::from(threads - 1);
        log::info!("Starting read recruitment with 1 read/write thread, and {} recruitment threads", n_workers);
        let mut main_worker = MainWorker::<T, _, _>::new(self, reader, writers, n_workers, self.params.chunk_size);
        main_worker.run()?;
        main_worker.finish()
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    pub fn recruit<T: RecruitableRecord, W: io::Write>(
        &self,
        reader: impl FastxRead<Record = T>,
        mut writers: Vec<W>,
        threads: u16,
    ) -> Result<(), Error>
    {
        assert_eq!(writers.len(), self.n_loci as usize, "Unexpected number of writers");
        if threads <= 1 {
            self.recruit_single_thread(reader, &mut writers)
        } else {
            self.recruit_multi_thread(reader, writers, threads)
        }
    }
}

/// Vector of records and corresponding answers.
/// Is send between threads: main thread reads records and sends the vector to workers.
/// In the meantime, workers receive records and fills corresponding answers (recruitment targets for the record),
/// and then send shipments back to the main thread.
type Shipment<T> = Vec<(T, Answer)>;

struct Worker<T> {
    targets: Targets,
    buffer1: Vec<Minimizer>,
    buffer2: IntMap<u16, u32>,
    /// Receives records that need to be recruited.
    receiver: Receiver<Shipment<T>>,
    /// Sends already recruited reads back to the main thread.
    sender: Sender<Shipment<T>>,
}

impl<T> Worker<T> {
    fn new(targets: Targets, receiver: Receiver<Shipment<T>>, sender: Sender<Shipment<T>>) -> Self {
        Self {
            buffer1: Vec::new(),
            buffer2: IntMap::default(),
            targets, receiver, sender,
        }
    }
}

impl<T: RecruitableRecord> Worker<T> {
    fn run(mut self) {
        // Block thread and wait for the shipment.
        while let Ok(mut shipment) = self.receiver.recv() {
            assert!(!shipment.is_empty());
            for (record, answer) in shipment.iter_mut() {
                record.recruit(&self.targets, answer, &mut self.buffer1, &mut self.buffer2);
            }
            if let Err(_) = self.sender.send(shipment) {
                log::error!("Read recruitment: main thread stopped before the child thread.");
                return;
            }
        }
    }
}

/// Worker in the main thread, that organizes other workers, as well as reads/writes reads.
struct MainWorker<T, R: FastxRead<Record = T>, W> {
    /// Fasta/q reader.
    /// Becomes `None`, once the stream has ended.
    reader: Option<R>,
    /// Fasta/q writers for each of the loci.
    writers: Vec<W>,
    /// Senders from the main thread to the workers. Sends reads to be analyzed.
    senders: Vec<Sender<Shipment<T>>>,
    /// Receivers from workers to the main thread. Receives
    /// - analyzed reads with possible recruited loci,
    /// - bool: true if any of the reads were recruited.
    receivers: Vec<Receiver<Shipment<T>>>,
    /// Thread handles.
    handles: Vec<thread::JoinHandle<()>>,

    /// Send reads between threads in vectors of this size.
    chunk_size: usize,
    /// Does the worker currently recruit reads?
    is_busy: Vec<bool>,
    /// Recruitment statistics.
    stats: Stats,

    /// Chunks of reads that were read from the reader and are ready to be analyzed.
    to_send: Vec<Shipment<T>>,
    /// Chunks of reads that were analyzed and need to be writter to the writers.
    to_write: Vec<Shipment<T>>,
}

impl<T, R, W> MainWorker<T, R, W>
where T: RecruitableRecord,
      R: FastxRead<Record = T>,
      W: io::Write,
{
    fn new(targets: &Targets, reader: R, writers: Vec<W>, n_workers: usize, chunk_size: usize) -> Self {
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        let mut handles = Vec::with_capacity(n_workers);

        for _ in 0..n_workers {
            let (sender1, receiver1) = mpsc::channel();
            let (sender2, receiver2) = mpsc::channel();
            let worker = Worker::new(targets.clone(), receiver1, sender2);
            senders.push(sender1);
            receivers.push(receiver2);
            handles.push(thread::spawn(|| worker.run()));
        }
        Self {
            writers, senders, receivers, handles, chunk_size,
            reader: Some(reader),
            is_busy: vec![false; n_workers],
            stats: Stats::new(),
            to_send: Vec::new(),
            to_write: Vec::new(),
        }
    }

    /// Starts the process: provides the first task to each worker.
    fn start(&mut self) -> Result<(), Error> {
        for (is_busy, sender) in self.is_busy.iter_mut().zip(&self.senders) {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size)?;
            if !shipment.is_empty() {
                *is_busy = true;
                sender.send(shipment).expect("Recruitment worker has failed!");
            }
            if self.reader.is_none() {
                break;
            }
        }
        Ok(())
    }

    fn recv_send_iteration(&mut self) -> bool {
        let mut any_action = false;
        for ((receiver, sender), is_busy) in self.receivers.iter().zip(&self.senders)
                .zip(self.is_busy.iter_mut()) {
            if *is_busy {
                match receiver.try_recv() {
                    Ok(recv_shipment) => {
                        any_action = true;
                        self.to_write.push(recv_shipment);
                        if let Some(send_shipment) = self.to_send.pop() {
                            sender.send(send_shipment).expect("Recruitment worker has failed!");
                        } else {
                            *is_busy = false;
                        }
                    }
                    Err(TryRecvError::Empty) => { continue; }
                    Err(TryRecvError::Disconnected) => panic!("Recruitment worker has failed!"),
                }
            } else if let Some(send_shipment) = self.to_send.pop() {
                any_action = true;
                sender.send(send_shipment).expect("Recruitment worker has failed!");
                *is_busy = true;
            }
        }
        any_action
    }

    fn write_read_iteration(&mut self) -> Result<(), Error> {
        if self.reader.is_none() {
            return Ok(())
        }
        while let Some(mut shipment) = self.to_write.pop() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
            self.stats.timed_print_log();
            fill_shipment(&mut self.reader, &mut shipment)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
            if self.reader.is_none() {
                return Ok(());
            }
        }
        if self.to_send.is_empty() {
            let shipment = read_new_shipment(&mut self.reader, self.chunk_size)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
        }
        Ok(())
    }

    /// Main part of the multi-thread recruitment.
    /// Iterates until there are any shipments left to read from the input files.
    fn run(&mut self) -> Result<(), Error> {
        self.start()?;
        while self.reader.is_some() || !self.to_send.is_empty() {
            self.write_read_iteration()?;
            // There were no updates, and there are shipments ready to be sent.
            if !self.recv_send_iteration() && !self.to_send.is_empty() {
                const SLEEP: Duration = Duration::from_micros(100);
                thread::sleep(SLEEP);
            }
        }
        Ok(())
    }

    /// Finish the main thread: write all remaining shipments to the output files, and stop worker threads.
    fn finish(mut self) -> Result<(), Error> {
        assert!(self.reader.is_none() && self.to_send.is_empty());
        for shipment in self.to_write.into_iter() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
        }
        for (&is_busy, receiver) in self.is_busy.iter().zip(&self.receivers) {
            if is_busy {
                // Block thread and wait for the task completion.
                let shipment = receiver.recv().expect("Recruitment worker has failed!");
                write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
                self.stats.timed_print_log();
            }
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason");
        }
        self.stats.finish();
        Ok(())
    }
}

/// Fills `shipment` from the reader.
/// Output shipment may be empty, if the stream has ended.
fn fill_shipment<T, R>(opt_reader: &mut Option<R>, shipment: &mut Shipment<T>) -> Result<(), Error>
where T: Default,
      R: FastxRead<Record = T>,
{
    let reader = opt_reader.as_mut().expect("fill_shipment: reader must not be None");
    let mut new_len = 0;
    for (record, _) in shipment.iter_mut() {
        if reader.read_next(record)? {
            new_len += 1;
        } else {
            shipment.truncate(new_len);
            *opt_reader = None;
            break;
        }
    }
    Ok(())
}

fn read_new_shipment<T, R>(opt_reader: &mut Option<R>, chunk_size: usize) -> Result<Shipment<T>, Error>
where T: Clone + Default,
      R: FastxRead<Record = T>,
{
    let mut shipment = vec![Default::default(); chunk_size];
    fill_shipment(opt_reader, &mut shipment)?;
    Ok(shipment)
}

/// Writes recruited records to the output files.
fn write_shipment<T>(writers: &mut [impl io::Write], shipment: &Shipment<T>, stats: &mut Stats) -> Result<(), Error>
where T: fastx::WritableRecord,
{
    stats.processed += shipment.len() as u64;
    for (record, answer) in shipment.iter() {
        stats.recruited += u64::from(!answer.is_empty());
        for &locus_ix in answer.iter() {
            record.write_to(&mut writers[usize::from(locus_ix)]).map_err(add_path!(!))?;
        }
    }
    Ok(())
}
