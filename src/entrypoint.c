// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_gpredomicsR_extendr(void *dll);

void R_init_gpredomicsR(void *dll) {
    R_init_gpredomicsR_extendr(dll);
}
