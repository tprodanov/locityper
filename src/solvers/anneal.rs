use std::fmt;
use rand::{
    Rng, SeedableRng,
    rngs::SmallRng,
};
use crate::{
    Error,
    model::assgn::ReadAssignment,
};
use super::Solver;
