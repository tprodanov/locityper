//! Serialize and deserialize data into JSON format.

use crate::Error;

/// The object can be serialized and deserialized into JSON.
pub trait JsonSer: Sized {
    fn save(&self) -> json::JsonValue;

    fn load(obj: &json::JsonValue) -> Result<Self, Error>;
}

pub fn parse_f64_arr(obj: &json::JsonValue, key: &str, arr: &mut [f64]) -> Result<(), Error> {
    if let json::JsonValue::Array(v) = &obj[key] {
        if v.len() != arr.len() {
            return Err(Error::JsonLoad(format!("Failed to parse '{}': incorrect number of elements in array '{}'",
                obj, key)));
        }
        for (i, val) in v.iter().enumerate() {
            arr[i] = val.as_f64().ok_or_else(|| Error::JsonLoad(
                format!("Failed to parse '{}': element #{} of array '{}' is not a float", obj, i, key)))?;
        }
        Ok(())
    } else {
        Err(Error::JsonLoad(format!("Failed to parse '{}': missing or incorrect array '{}'", obj, key)))
    }
}

#[doc(hidden)]
pub fn json_err(obj: &json::JsonValue, key: &str) -> Error {
    let mut obj_str = obj.to_string();
    const SUBSTR_SIZE: usize = 200;
    if obj_str.len() > SUBSTR_SIZE {
        obj_str = format!("{}...", &obj_str[..SUBSTR_SIZE]);
    }
    return Error::JsonLoad(
        format!("Failed to parse '{}': missing or incorrect '{}' field!", obj_str, key));
}

/// Function that simplifies fetching from Json value.
/// Usage:
/// ```
/// // Creates three variables and fetches their values (var1: f64, var2: usize and var3: bool).
/// json_get!(obj => var1 (as_f64), var2 (as_usize), var3 (as_bool));
/// ```
/// The macro also supports optional values, using `json_get!(obj => var1? (as_f64))`
/// will set `var1` to `Option<f64>`.
/// Nevertheless, the error is produced if the variable is present in the json file, but has incorrect type.
macro_rules! json_get {
    ($obj:expr => $var:ident ($convert:ident)) => {
        let $var = match $obj[std::stringify!($var)].$convert() {
            Some(val) => val,
            None => Err($crate::bg::ser::json_err($obj, stringify!($var)))?,
        };
    };
    ($obj:expr => $var:ident ? ($convert:ident)) => {
        let val = &$obj[std::stringify!($var)];
        let $var = if val.is_null() {
            None
        } else {
            match $obj[std::stringify!($var)].$convert() {
                Some(val) => Some(val),
                None => Err($crate::bg::ser::json_err($obj, stringify!($var)))?,
            }
        };
    };
    ($obj:expr => $var:ident ($convert:ident) , $($tail:tt)*) => {
        json_get!($obj => $var ($convert));
        json_get!($obj => $($tail)*);
    };
    ($obj:expr => $var:ident ? ($convert:ident) , $($tail:tt)*) => {
        json_get!($obj => $var ? ($convert));
        json_get!($obj => $($tail)*);
    };
}
pub(crate) use json_get;
