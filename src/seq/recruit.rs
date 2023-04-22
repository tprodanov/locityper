//! Recruit reads from FASTQ files according to their k-mers.

use std::{
    io, thread,
    collections::hash_map::Entry,
    time::Instant,
    sync::mpsc::{self, Sender, Receiver, TryRecvError},
};
use const_format::formatcp;
use nohash::IntMap;
use smallvec::{smallvec, SmallVec};
use crate::{
    seq::{
        ContigSet,
        kmers,
        fastx::{RecordExt, FastxRead},
    },
};

pub struct Params {
    /// Recruit reads using k-mers of size `minimizer_k` that have the minimial hash across `minimizer_w`
    /// consecutive k-mers.
    /// Default: 15 and 10.
    pub minimizer_k: u8,
    pub minimizer_w: u8,
    /// Recruit reads that have at least this number of matching minimizers with one of the targets. Default: 2.
    pub min_matches: u8,
    /// Pass reads between threads in chunks of this size. Default: 10000.
    pub chunk_size: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            minimizer_k: 15,
            minimizer_w: 10,
            min_matches: 2,
            chunk_size: 10000,
        }
    }
}

/// Recruitment statistics: how long did it take, how many reads processed and how many recruited.
struct Stats {
    timer: Instant,
    recruited: u64,
    processed: u64,
    /// Last log message was at this number of processed reads.
    last_processed: u64,
}

impl Stats {
    fn new() -> Self {
        Self {
            timer: Instant::now(),
            recruited: 0,
            processed: 0,
            last_processed: 0,
        }
    }

    /// Prints log message if the number of processed reads divides `modulo`.
    fn print_log_mod(&mut self, modulo: u64) {
        if self.processed % modulo == 0 {
            self.print_log();
        }
    }

    /// Prints log message if there are more than `n_new` processed reads.
    fn print_log_diff(&mut self, n_new: u64) {
        if self.processed - self.last_processed >= n_new {
            self.print_log();
        }
    }

    fn print_log(&mut self) {
        let per_read = self.timer.elapsed().as_secs_f64() / self.processed as f64;
        log::debug!("    Recruited {:11} /{:11} reads,  {:.2} us/read", self.recruited, self.processed,
            per_read * 1e6);
        self.last_processed = self.processed;
    }

    fn finish(&mut self) {
        if self.processed > self.last_processed {
            self.print_log();
        }
        log::info!("Finished recruitment in {}", crate::ext::fmt::Duration(self.timer.elapsed()));
    }
}

const CAPACITY: usize = 4;
/// Key: minimizer, value: vector of loci indices, where the minimizer appears.
type MinimToLoci = IntMap<u32, SmallVec<[u16; CAPACITY]>>;
/// Vector of loci indices, to which the read was recruited.
type Answer = Vec<u16>;

/// Recruitement targets.
#[derive(Clone)]
pub struct Targets {
    minimizer_k: u8,
    minimizer_w: u8,
    /// Minimizers appearing across the targets.
    minim_to_loci: MinimToLoci,
    /// Minimal number of k-mer matches per read/read-pair, needed to recruit the read to the corresponding locus.
    min_matches: u8,
    n_loci: u16,
}

impl Targets {
    pub fn new<'a>(sets: impl Iterator<Item = &'a ContigSet>, params: &Params) -> Self {
        log::info!("Generating recruitment targets");
        let mut minim_to_loci = MinimToLoci::default();
        let mut buffer = Vec::new();

        let mut n_seqs = 0;
        let mut locus_ix: u16 = 0;
        for contig_set in sets {
            for contig_id in contig_set.contigs().ids() {
                buffer.clear();
                kmers::minimizers(contig_set.get_seq(contig_id), params.minimizer_k, params.minimizer_w, &mut buffer);
                for &minimizer in buffer.iter() {
                    match minim_to_loci.entry(minimizer) {
                        Entry::Occupied(entry) => {
                            let loci_ixs = entry.into_mut();
                            if *loci_ixs.last().unwrap() != locus_ix {
                                loci_ixs.push(locus_ix);
                            }
                        }
                        Entry::Vacant(entry) => { entry.insert(smallvec![locus_ix]); }
                    }
                }
            }
            locus_ix = locus_ix.checked_add(1)
                .expect(formatcp!("Too many contig sets (allowed at most {})", u16::MAX));
            n_seqs += contig_set.len();
        }
        let n_loci = locus_ix;
        log::info!("Collected {} minimizers across {} loci and {} sequences", minim_to_loci.len(), n_loci, n_seqs);
        Self {
            minim_to_loci, n_loci,
            minimizer_k: params.minimizer_k,
            minimizer_w: params.minimizer_w,
            min_matches: params.min_matches,
        }
    }

    /// Record a specific single- or paired-read to one or more loci.
    ///
    /// The read is written to all loci (to the corresponding `out_files[locus_ix]`),
    /// for which the number of loci_matches will be at least `min_matches`.
    /// Returns true if the read was recruited anywhere.
    fn recruit_record<T: RecordExt>(
        &self,
        record: &T,
        answer: &mut Answer,
        buffer_minimizers: &mut Vec<u32>,
        buffer_matches: &mut IntMap<u16, u8>,
    ) {
        buffer_minimizers.clear();
        record.minimizers(self.minimizer_k, self.minimizer_w, buffer_minimizers);
        buffer_matches.clear();
        for minimizer in buffer_minimizers.iter() {
            if let Some(loci_ixs) = self.minim_to_loci.get(minimizer) {
                for &locus_ix in loci_ixs.iter() {
                    buffer_matches.entry(locus_ix)
                        .and_modify(|counter| *counter = counter.saturating_add(1))
                        .or_insert(1);
                }
            }
        }
        answer.clear();
        for (&locus_ix, &count) in buffer_matches.iter() {
            if count >= self.min_matches {
                answer.push(locus_ix);
            }
        }
    }

    fn recruit_single_thread<T: RecordExt, W: io::Write>(
        &self,
        mut reader: impl FastxRead<Record = T>,
        writers: &mut [W],
    ) -> io::Result<()>
    {
        let mut record = T::default();
        let mut answer = Answer::new();
        let mut buffer1 = Vec::new();
        let mut buffer2 = IntMap::default();

        let mut stats = Stats::new();
        while reader.read_next(&mut record)? {
            self.recruit_record(&record, &mut answer, &mut buffer1, &mut buffer2);
            for &locus_ix in answer.iter() {
                record.write_to(&mut writers[usize::from(locus_ix)])?;
            }
            stats.recruited += u64::from(!answer.is_empty());
            stats.processed += 1;
            stats.print_log_mod(100_000);
        }
        stats.finish();
        Ok(())
    }

    fn recruit_multi_thread<T: RecordExt, W: io::Write>(
        &self,
        reader: impl FastxRead<Record = T>,
        writers: Vec<W>,
        threads: u16,
        chunk_size: usize,
    ) -> io::Result<()> {
        let n_workers = usize::from(threads - 1);
        log::info!("Starting read recruitment with 1 read/write thread, and {} recruitment threads", n_workers);
        let mut main_worker = MainWorker::<T, _, _>::new(self, reader, writers, n_workers, chunk_size);
        main_worker.run()?;
        main_worker.finish()
    }

    /// Recruit reads to the targets, possibly in multiple threads.
    pub fn recruit<T: RecordExt, W: io::Write>(
        &self,
        reader: impl FastxRead<Record = T>,
        mut writers: Vec<W>,
        threads: u16,
        chunk_size: usize,
    ) -> io::Result<()>
    {
        assert_eq!(writers.len(), self.n_loci as usize, "Unexpected number of writers");
        if threads <= 2 {
            if threads == 2 {
                log::warn!("2-thread recruitment is slower, running in single-thread mode!");
            }
            self.recruit_single_thread(reader, &mut writers)
        } else {
            self.recruit_multi_thread(reader, writers, threads, chunk_size)
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
    buffer1: Vec<u32>,
    buffer2: IntMap<u16, u8>,
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

impl<T: RecordExt> Worker<T> {
    fn run(mut self) {
        // Block thread and wait for the shipment.
        while let Ok(mut shipment) = self.receiver.recv() {
            assert!(!shipment.is_empty());
            for (record, answer) in shipment.iter_mut() {
                self.targets.recruit_record(record, answer, &mut self.buffer1, &mut self.buffer2);
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
    in_progress: Vec<bool>,
    /// Recruitment statistics.
    stats: Stats,

    /// Chunks of reads that were read from the reader and are ready to be analyzed.
    to_send: Vec<Shipment<T>>,
    /// Chunks of reads that were analyzed and need to be writter to the writers.
    to_write: Vec<Shipment<T>>,
}

impl<T, R, W> MainWorker<T, R, W>
where T: RecordExt,
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
            handles.push(thread::spawn(move || worker.run()));
        }
        Self {
            writers, senders, receivers, handles, chunk_size,
            reader: Some(reader),
            in_progress: vec![false; n_workers],
            stats: Stats::new(),
            to_send: Vec::new(),
            to_write: Vec::new(),
        }
    }

    /// Starts the process: provides the first task to each worker.
    fn start(&mut self) -> io::Result<()> {
        for (in_progress, sender) in self.in_progress.iter_mut().zip(&self.senders) {
            let mut shipment = vec![Default::default(); self.chunk_size];
            fill_shipment(&mut self.reader, &mut shipment)?;
            if !shipment.is_empty() {
                *in_progress = true;
                sender.send(shipment).map_err(|_| io::Error::new(
                    io::ErrorKind::Other, "Recruitment worker has failed!"))?;
            }
            if self.reader.is_none() {
                break;
            }
        }
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        assert!(self.reader.is_none() && self.to_send.is_empty());
        for shipment in self.to_write.into_iter() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
        }
        for (&in_progress, receiver) in self.in_progress.iter().zip(&self.receivers) {
            if in_progress {
                // Block thread and wait for the task completion.
                let shipment = receiver.recv().map_err(|_| io::Error::new(
                    io::ErrorKind::Other, "Recruitment worker has failed!"))?;
                write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
            }
        }
        std::mem::drop(self.senders);
        for handle in self.handles.into_iter() {
            handle.join().expect("Process failed for unknown reason");
        }
        self.stats.finish();
        Ok(())
    }

    fn try_recv_iteration(&mut self) -> io::Result<()> {
        for (in_progress, receiver) in self.in_progress.iter_mut().zip(&self.receivers) {
            if *in_progress {
                match receiver.try_recv() {
                    Ok(shipment) => {
                        *in_progress = false;
                        self.to_write.push(shipment);
                    }
                    Err(TryRecvError::Empty) => {}
                    Err(TryRecvError::Disconnected) => return Err(io::Error::new(
                        io::ErrorKind::Other, "Recruitment worker has failed!")),
                }
            }
        }
        Ok(())
    }

    fn write_read_iteration(&mut self) -> io::Result<()> {
        if self.reader.is_none() {
            return Ok(())
        }
        while let Some(mut shipment) = self.to_write.pop() {
            write_shipment(&mut self.writers, &shipment, &mut self.stats)?;
            fill_shipment(&mut self.reader, &mut shipment)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
            if self.reader.is_none() {
                return Ok(());
            }
        }
        if self.to_send.is_empty() {
            let mut shipment = vec![Default::default(); self.chunk_size];
            fill_shipment(&mut self.reader, &mut shipment)?;
            if !shipment.is_empty() {
                self.to_send.push(shipment);
            }
        }
        Ok(())
    }

    fn send_iteration(&mut self) -> io::Result<()> {
        for (in_progress, sender) in self.in_progress.iter_mut().zip(&self.senders) {
            if !*in_progress {
                if let Some(shipment) = self.to_send.pop() {
                    *in_progress = true;
                    sender.send(shipment).map_err(|_| io::Error::new(
                        io::ErrorKind::Other, "Recruitment worker has failed!"))?;
                } else {
                    break;
                }
            }
        }
        Ok(())
    }

    fn run(&mut self) -> io::Result<()> {
        self.start()?;
        while self.reader.is_some() || !self.to_send.is_empty() {
            self.try_recv_iteration()?;
            self.write_read_iteration()?;
            self.send_iteration()?;
            self.stats.print_log_diff(1_000_000);
        }
        Ok(())
    }
}

/// Fills `shipment` from the reader.
/// Output shipment may be empty, if the stream has ended.
fn fill_shipment<T, R>(opt_reader: &mut Option<R>, shipment: &mut Shipment<T>) -> io::Result<()>
where T: RecordExt,
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

/// Writes recruited records to the output files.
fn write_shipment<T, W>(writers: &mut [W], shipment: &Shipment<T>, stats: &mut Stats) -> io::Result<()>
where T: RecordExt,
      W: io::Write,
{
    stats.processed += shipment.len() as u64;
    for (record, answer) in shipment.iter() {
        stats.recruited += u64::from(!answer.is_empty());
        for &locus_ix in answer.iter() {
            record.write_to(&mut writers[usize::from(locus_ix)])?;
        }
    }
    Ok(())
}
