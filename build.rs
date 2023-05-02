use std::{
    env,
    path::PathBuf,
};

/// WFA2-lib bindings.
fn wfa() {
    const WFA_PATH: &'static str = "WFA2";

    // The directory of the WFA libraries, added to the search path.
    println!("cargo:rustc-link-search={}/lib", WFA_PATH);
    // Link the `wfa-lib` library.
    println!("cargo:rustc-link-lib=wfa");
    // Invalidate the built crate whenever the linked library changes.
    println!("cargo:rerun-if-changed={}/lib/libwfa.a", WFA_PATH);

    // Generate bindings.
    let bindings = bindgen::Builder::default()
        // Generate bindings for this header file.
        .header(&format!("{}/utils/commons.h", WFA_PATH))
        .header(&format!("{}/wavefront/wavefront_align.h", WFA_PATH))
        .header(&format!("{}/alignment/cigar.h", WFA_PATH))
        // Add this directory to the include path to find included header files.
        .clang_arg(&format!("-I{}", WFA_PATH))
        // Generate bindings for all functions starting with `wavefront_`.
        .allowlist_function("wavefront_.*")
        // .allowlist_function("cigar_sprint")
        // Generate bindings for all variables starting with `wavefront_`.
        .allowlist_var("wavefront_.*")
        // Invalidate the built crate whenever any of the included header files
        // changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings_wfa.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings_wfa.rs"))
        .expect("Couldn't write bindings!");
}

fn main() {
    wfa();
}
