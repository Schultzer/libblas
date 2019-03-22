use libblas::level3;
use num_complex::Complex;
use rand::Rng;
use test::Bencher;

#[bench]
fn gemm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::gemm(
            'n',
            'n',
            6,
            6,
            4,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m.clone(),
            6,
            Complex::new(0.2, 0.8),
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn hemm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::hemm(
            'l',
            'u',
            6,
            6,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m.clone(),
            6,
            Complex::new(0.2, 0.8),
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn symm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::symm(
            'l',
            'u',
            6,
            6,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m.clone(),
            6,
            Complex::new(0.2, 0.8),
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn syrk(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::syrk(
            'l',
            't',
            6,
            3,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            Complex::new(0.2, 0.8),
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn syr2k(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::syr2k(
            'l',
            't',
            6,
            3,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m.clone(),
            6,
            Complex::new(0.2, 0.8),
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn herk(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::herk(
            'u',
            'n',
            6,
            3,
            0.8,
            &m.clone(),
            6,
            0.2,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn her2k(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::her2k(
            'l',
            'n',
            6,
            3,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m.clone(),
            6,
            0.2,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn trmm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::trmm(
            'l',
            'u',
            'n',
            'u',
            3,
            0,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn trsm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level3::complex::trsm(
            'l',
            'u',
            'n',
            'u',
            3,
            0,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &mut m.clone(),
            6,
        ))
    })
}
