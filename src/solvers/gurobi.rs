use std::{fmt, result::Result};
use grb::{
    *,
    expr::{LinExpr, GurobiSum},
};
use crate::{
    Error,
    model::{
        windows::UNMAPPED_WINDOW,
        assgn::ReadAssignment,
    },
};
use super::Solver;
