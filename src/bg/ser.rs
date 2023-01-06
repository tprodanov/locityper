//! Serialize and deserialize data into JSON format.

use std::fmt::{self, Display, Formatter};

/// Loading error.
#[derive(Clone, Debug)]
pub struct LoadError(pub String);

impl Display for LoadError {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// The object can be serialized and deserialized into JSON.
pub trait JsonSer: Sized {
    fn save(&self) -> json::JsonValue;

    fn load(obj: &json::JsonValue) -> Result<Self, LoadError>;
}
