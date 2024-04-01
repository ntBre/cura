fn main() {
    let out = std::env::var("DEP_SHIM2_INCLUDE").unwrap();
    println!("cargo:rustc-link-arg=-Wl,-rpath,{out}",);
}
