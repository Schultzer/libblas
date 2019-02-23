#[macro_use]
extern crate criterion;
extern crate blas;
extern crate blas_src;

use blas::*;
use criterion::Criterion;
mod fixtures;

fn blas_src_gemm() -> () {
    let mat = fixtures::matrix_mxn(4, 4);
    let mut c = mat.clone();
    unsafe {
        dgemm(b'n', b'n', 4, 4, 4, 0.3, &mat, 4, &mat, 4, -1.2, &mut c, 4);
    }
}

fn blas_src_axpy() -> () {
    let x = vec![
        1.2629542848807933098,
        -0.3262333607056494000,
        1.3297992629225006134,
        1.2724293214294046805,
        0.4146414344564082199,
        -1.5399500419037095433,
        -0.9285670347135380753,
        -0.2947204467905601977,
        -0.0057671727475369552,
        2.4046533888579508798,
    ];
    let mut y = vec![
        -0.560475646552212603,
        -0.230177489483279957,
        1.558708314149124030,
        0.070508391424576003,
        0.129287735160946243,
        1.715064986883281017,
        0.460916205989202299,
        -1.265061234606533969,
        -0.686852851893526073,
        -0.445661970099958060,
    ];
    unsafe {
        daxpy(10, 23.0, &x, 1, &mut y, 1);
    }
}

fn axpy() -> () {
    let x = vec![
        1.2629542848807933098,
        -0.3262333607056494000,
        1.3297992629225006134,
        1.2724293214294046805,
        0.4146414344564082199,
        -1.5399500419037095433,
        -0.9285670347135380753,
        -0.2947204467905601977,
        -0.0057671727475369552,
        2.4046533888579508798,
    ];
    let mut y = vec![
        -0.560475646552212603,
        -0.230177489483279957,
        1.558708314149124030,
        0.070508391424576003,
        0.129287735160946243,
        1.715064986883281017,
        0.460916205989202299,
        -1.265061234606533969,
        -0.686852851893526073,
        -0.445661970099958060,
    ];
    blasrs::level1::axpy(10, 23.0, &x, 1, &mut y, 1);
}

fn gemm() -> () {
    let mat = fixtures::matrix_mxn(6, 6);
    blasrs::level3::gemm(
        'n',
        'n',
        4,
        6,
        3,
        0.3,
        &mat,
        6,
        &mat,
        6,
        -1.2,
        &mut mat.clone(),
        6,
    );
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("blasrs::level3::gemm", |b| b.iter(|| gemm()));
    c.bench_function("blasrs::level1::axpy", |b| b.iter(|| axpy()));
    c.bench_function("blas-src accelerate dgemm", |b| b.iter(|| blas_src_gemm()));
    c.bench_function("blas-src accelerate daxpy", |b| b.iter(|| blas_src_axpy()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
