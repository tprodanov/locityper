use std::{
    io, process,
};

/// RAII child wrapper, that kills the child if it gets dropped.
pub struct ChildGuard {
    child: process::Child,
    armed: bool,
}

impl ChildGuard {
    pub fn new(child: process::Child) -> Self {
        Self {
            child,
            armed: true,
        }
    }

    pub fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for ChildGuard {
    fn drop(&mut self) {
        if self.armed {
            match self.child.kill() {
                Err(e) => {
                    // InvalidInput means that the process exited already.
                    if e.kind() != io::ErrorKind::InvalidInput {
                        log::error!("Could not kill child process: {}", e);
                    }
                }
                Ok(_) => log::error!("Successfully killed child process"),
            }
        }
    }
}