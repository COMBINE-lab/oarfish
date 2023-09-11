fn main() {
    println!("cargo:rustc-link-search=/mnt/scratch4/zahra/oarfish/src");
    println!("cargo:rustc-link-lib=mlpack_kde_wrapper");
}