use std::{
    path::PathBuf,
    sync::Arc,
    io::{self, Write},
    ops::{Deref, DerefMut},
};
use crate::{
    ext::{
        self,
        sys::CreatableFile,
    },
    err::{error, add_path},
};

/// A writer that keeps the files alive until it is dropped.
pub struct ThreadFile<W> {
    file: W,
    _guard: Arc<()>,
}

impl<W> ThreadFile<W> {
    #[inline]
    fn new(file: W, guard: Arc<()>) -> Self {
        Self {
            file,
            _guard: guard,
        }
    }
}

impl<W> Deref for ThreadFile<W> {
    type Target = W;

    #[inline(always)]
    fn deref(&self) -> &W { &self.file }
}

impl<W> DerefMut for ThreadFile<W> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut W { &mut self.file }
}

impl<W: Write> Write for ThreadFile<W> {
    #[inline(always)]
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> { self.file.write(buf) }

    #[inline(always)]
    fn flush(&mut self) -> io::Result<()> { self.file.flush() }
}

pub struct ThreadFiles {
    /// Final filename. Will also be used by the first thread.
    out_filename: Option<PathBuf>,
    /// Temporary filenames in multi-threaded application.
    tmp_filenames: Vec<PathBuf>,
    /// File under this index will be created next.
    curr_index: usize,
    /// If present, write this text into the first file.
    header: Option<String>,
    threads: usize,
    guard: Arc<()>,
}

impl ThreadFiles {
    pub fn new(out_filename: PathBuf, threads: usize, header: Option<String>) -> Self {
        let extension = out_filename.extension();
        let tmp_filenames = (2..=threads).map(|i| {
            let mut tmp_filename = out_filename.clone();
            tmp_filename.set_extension(format!("tmp{}", i));
            if let Some(e) = extension {
                tmp_filename.add_extension(e);
            }
            tmp_filename
        }).collect();
        Self::with_filenames(out_filename, tmp_filenames, header)
    }

    /// Creates thread files with the first filename and temporary filenames for thread 2, 3, etc.
    #[inline]
    pub fn with_filenames(
        out_filename: PathBuf,
        tmp_filenames: Vec<PathBuf>,
        header: Option<String>,
    ) -> Self {
        Self {
            threads: 1 + tmp_filenames.len(),
            tmp_filenames, header,
            out_filename: Some(out_filename),
            curr_index: 0,
            guard: Arc::new(()),
        }
    }

    /// Creates disabled thread files (no files will be created).
    pub fn disabled(threads: usize) -> Self {
        Self {
            out_filename: None,
            tmp_filenames: Vec::new(),
            curr_index: 0,
            header: None,
            threads,
            guard: Arc::new(()),
        }
    }

    /// Creates next file. Returns None if files were disabled or too many were already constructed.
    pub fn create_optional_next<W: CreatableFile>(&mut self) -> crate::Result<Option<ThreadFile<W>>> {
        if self.curr_index >= self.len() {
            return Ok(None);
        }
        let file = if self.curr_index == 0 {
            let filename = self.out_filename.as_ref().expect("Temporary files must not be disabled here");
            let mut file = W::create(filename)?;
            if let Some(text) = &self.header {
                writeln!(file, "{}", text).map_err(add_path!(filename))?;
            }
            file
        } else {
            let filename = &self.tmp_filenames[self.curr_index - 1];
            W::create(filename)?
        };
        self.curr_index += 1;
        Ok(Some(ThreadFile::new(file, Arc::clone(&self.guard))))
    }

    // /// Same as `create_optional_next`, but panics if there is no next file.
    // pub fn create_next<W: CreatableFile>(&mut self) -> crate::Result<ThreadFile<W>> {
    //     let i = self.curr_index;
    //     self.create_optional_next()?
    //         .ok_or_else(|| error!(RuntimeError, "Requested {}-th temporary file ({} available)", i + 1, self.len()))
    // }

    /// Creates Option<File> in an iterator, will return None's if these thread files were disabled.
    pub fn iter_options<W: CreatableFile>(&mut self) -> impl Iterator<Item = crate::Result<Option<ThreadFile<W>>>> {
        (self.curr_index..self.threads).map(|_| self.create_optional_next())
    }

    // /// Output next files as an iterator.
    // pub fn iter<W: CreatableFile>(&mut self) -> impl Iterator<Item = crate::Result<ThreadFile<W>>> {
    //     (self.curr_index..self.threads).map(|_| self.create_next())
    // }

    /// Merge temporary files into the first one, and delete them.
    pub fn merge(&mut self) -> crate::Result<()> {
        if let Some(out_filename) = &self.out_filename {
            if Arc::strong_count(&self.guard) > 1 {
                return Err(error!(RuntimeError,
                    "Tried to merge temporary files ({:?}) while some writers are still open", self.out_filename));
            }
            ext::sys::merge_files(out_filename, &self.tmp_filenames)?;
            self.tmp_filenames.clear();
        } else {
            debug_assert!(self.tmp_filenames.is_empty());
        }
        Ok(())
    }

    pub fn len(&self) -> usize {
        usize::from(self.out_filename.is_some()) + self.tmp_filenames.len()
    }
}

impl Drop for ThreadFiles {
    fn drop(&mut self) {
        if let Err(e) = self.merge() {
            log::error!("Error `{}` when dropping temporary files", e.display());
        }
    }
}
