use std::{
    io::{self, Read},
};
use brotli::{
    BrotliState, BrotliResult,
    BrotliDecompressStream as decompress_stream,
    enc::StandardAlloc,
};

use super::BUFFER_SIZE_4MB;

#[inline]
fn new_state() -> BrotliState<StandardAlloc, StandardAlloc, StandardAlloc> {
    BrotliState::new(Default::default(), Default::default(), Default::default())
}

/// Wrapper around brotli Decompressor that supports multiple concatenated streams.
pub struct MultiBrDecoder<R: Read> {
    stream: R,
    buf: Vec<u8>,
    /// Current buffer offset.
    offset: usize,
    /// How many bytes of the buffer were filled.
    rem_buf: usize,
    state: Option<BrotliState<StandardAlloc, StandardAlloc, StandardAlloc>>,
    last_result: BrotliResult,
}

impl<R: Read> MultiBrDecoder<R> {
    pub fn new(stream: R) -> Self {
        Self {
            stream,
            buf: vec![0; BUFFER_SIZE_4MB],
            offset: 0,
            rem_buf: 0,
            state: None,
            last_result: BrotliResult::ResultSuccess,
        }
    }

    #[inline]
    fn fill_buf(&mut self) -> io::Result<bool> {
        if self.rem_buf == 0 {
            self.offset = 0;
            self.rem_buf = self.stream.read(&mut self.buf)?;
            Ok(self.rem_buf > 0)
        } else {
            Ok(true)
        }
    }
}

impl<R: Read> Read for MultiBrDecoder<R> {
    fn read(&mut self, out_buf: &mut [u8]) -> io::Result<usize> {
        let mut out_offset = 0;
        let mut rem_out = out_buf.len();
        loop {
            // Have to check this way because BrotliResult implements neither Clone nor PartialEq
            if let &BrotliResult::NeedsMoreOutput = &self.last_result { /* Do nothing */ }
            else if !self.fill_buf()? {
                match self.last_result {
                    BrotliResult::ResultSuccess => return Ok(0),
                    BrotliResult::NeedsMoreInput => return Err(io::Error::new(io::ErrorKind::UnexpectedEof,
                        "truncated brotli stream")),
                    _ => unreachable!(),
                }
            }

            self.last_result = decompress_stream(
                &mut self.rem_buf, &mut self.offset, &self.buf,
                &mut rem_out, &mut out_offset, out_buf, &mut 0,
                self.state.get_or_insert_with(new_state));

            match self.last_result {
                BrotliResult::ResultSuccess => self.state = None,
                BrotliResult::NeedsMoreOutput => {}
                BrotliResult::NeedsMoreInput => continue,
                BrotliResult::ResultFailure =>
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "corrupt brotli stream")),
            }
            if out_offset > 0 {
                return Ok(out_offset);
            }
        }
    }
}
