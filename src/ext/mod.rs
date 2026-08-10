//! Extending standard functions and structures.

pub mod vec;
pub mod trimat;
pub mod fmt;
pub mod sys;
pub mod rand;

pub use trimat::TriangleMatrix;

enum State<T, F> {
    Uninit(F),
    Init(T),
    Err,
}

/// Lazy storage similar to LazyCell, but with a function that can fail.
/// On fail, error is returned once and subsequent requests will panic.
pub struct LazyResult<T, E, F: FnOnce() -> Result<T, E>> {
    state: State<T, F>,
}

impl<T, E, F: FnOnce() -> Result<T, E>> LazyResult<T, E, F> {
    pub fn new(constructor: F) -> Self {
        Self {
            state: State::Uninit(constructor),
        }
    }

    pub fn get(&mut self) -> Result<&T, E> {
        if let State::Uninit(_) = self.state {
            let State::Uninit(f) = std::mem::replace(&mut self.state, State::Err) else { unreachable!() };
            self.state = State::Init(f()?);
        }

        match &self.state {
            State::Init(val) => Ok(val),
            State::Uninit(_) => unreachable!(),
            State::Err => panic!("LazyResult was requested after an error"),
        }
    }
}
