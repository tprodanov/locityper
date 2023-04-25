use std::fmt;
use highs::{RowProblem, Model, Col, Sense, HighsModelStatus as Status};
use crate::{
    Error,
    model::{
        windows::UNMAPPED_WINDOW,
        assgn::ReadAssignment,
    },
    ext::vec::F64Ext,
};
use super::Solver;
