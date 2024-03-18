#[cfg(feature = "align")]
mod build_wfa {
    use std::{
        env,
        path::PathBuf,
    };

    const WFA_PATH: &'static str = "WFA2";

    fn generate_wfa_bindings() -> Result<bindgen::Bindings, bindgen::BindgenError> {
        bindgen::Builder::default()
            // Generate bindings for this header file.
            .header(&format!("{}/utils/commons.h", WFA_PATH))
            .header(&format!("{}/wavefront/wavefront_align.h", WFA_PATH))
            .header(&format!("{}/alignment/cigar.h", WFA_PATH))
            // Add this directory to the include path to find included header files.
            .clang_arg(&format!("-I{}", WFA_PATH))

            // Generate bindings for all functions and variables starting with `wavefront_`.
            .allowlist_function("wavefront_.*")
            .allowlist_var("wavefront_.*")
            // Invalidate the built crate whenever any of the included header files changed.
            .parse_callbacks(Box::new(bindgen::CargoCallbacks))
            // Finish the builder and generate the bindings.
            .generate()
    }

    fn make_wfa() {
        let output = std::process::Command::new("make")
            .arg("clean")
            .arg("all")
            .arg("BUILD_WFA_PARALLEL=0") // Disable parallelization.
            .args(&["CC_FLAGS+=-Wall", "CC_FLAGS+=-g", "CC_FLAGS+=-fPIC"])
            .current_dir(&WFA_PATH)
            .output()
            .expect("Failed to compile WFA2-lib");
        if !output.status.success() {
            panic!("Failed to compile WFA2-libs: {}", String::from_utf8_lossy(&output.stderr));
        }
    }

    /// WFA2-lib bindings.
    pub fn build() {
        // The directory of the WFA libraries, added to the search path.
        println!("cargo:rustc-link-search={}/lib", WFA_PATH);
        // Link the `wfa-lib` library.
        println!("cargo:rustc-link-lib=wfa");
        // Invalidate the built crate whenever the linked library changes.
        println!("cargo:rerun-if-changed={}/lib/libwfa.a", WFA_PATH);

        if !PathBuf::from(WFA_PATH).join("build/wavefront.o").exists() {
            make_wfa();
        }
        let bindings = generate_wfa_bindings().expect("Cannot generate WFA2-lib bindings");

        // Write the bindings to the $OUT_DIR/bindings_wfa.rs file.
        let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
        bindings.write_to_file(out_path.join("bindings_wfa.rs"))
            .expect("Couldn't write bindings!");
    }
}

fn main() {
    #[cfg(feature = "align")]
    { build_wfa::build(); }
}
