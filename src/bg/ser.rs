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

pub fn parse_f64_arr(obj: &json::JsonValue, key: &str, arr: &mut [f64]) -> Result<(), LoadError> {
    if let json::JsonValue::Array(v) = &obj[key] {
        if v.len() != arr.len() {
            return Err(LoadError(format!("Failed to parse '{}': incorrect number of elements in array '{}'",
                obj, key)));
        }
        for (i, val) in v.iter().enumerate() {
            arr[i] = val.as_f64().ok_or_else(||
                LoadError(format!("Failed to parse '{}': element #{} of array '{}' is not a float", obj, i, key)))?;
        }
        Ok(())
    } else {
        Err(LoadError(format!("Failed to parse '{}': missing or incorrect array '{}'", obj, key)))
    }
}
