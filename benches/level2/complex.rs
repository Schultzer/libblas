use libblas::level2;
use num_complex::Complex;
use rand::Rng;
use test::Bencher;

#[bench]
fn gbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::gbmv(
            'n',
            6,
            6,
            1,
            1,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m2.clone(),
            1,
            Complex::new(0.2, 0.8),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn gemv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::gemv(
            'n',
            6,
            6,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m2.clone(),
            1,
            Complex::new(0.2, 0.8),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn gerc(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::gerc(
            6,
            6,
            Complex::new(0.2, -0.8),
            &m2.clone(),
            1,
            &m2.clone(),
            1,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn geru(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::geru(
            6,
            6,
            Complex::new(0.2, -0.8),
            &m2.clone(),
            1,
            &m2.clone(),
            1,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn hbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::hbmv(
            'u',
            6,
            3,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m2.clone(),
            1,
            Complex::new(2.0, 5.0),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn hemv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::hemv(
            'u',
            6,
            Complex::new(0.2, 0.8),
            &m.clone(),
            6,
            &m2.clone(),
            1,
            Complex::new(2.0, 5.0),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn her(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }

    bh.iter(|| {
        test::black_box(level2::complex::her(
            'u',
            6,
            1.0,
            &m2.clone(),
            1,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn her2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::her2(
            'u',
            6,
            Complex::new(0.2, 0.8),
            &m2.clone(),
            1,
            &m2.clone(),
            1,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn hpmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::hpmv(
            'u',
            6,
            Complex::new(0.2, -0.8),
            &m,
            &m2.clone(),
            1,
            Complex::new(0.3, -0.7),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn hpr(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..64 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..8 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::hpr(
            'u',
            8,
            1.0,
            &m2.clone(),
            1,
            &mut m.clone(),
        ))
    })
}

#[bench]
fn hpr2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..64 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..8 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::hpr2(
            'u',
            8,
            Complex::new(0.2, -0.8),
            &m2.clone(),
            1,
            &m2.clone(),
            1,
            &mut m.clone(),
        ))
    })
}

#[bench]
fn tbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::tbmv(
            'u',
            'n',
            'n',
            6,
            5,
            &m.clone(),
            6,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn tbsv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::tbsv(
            'u',
            'n',
            'n',
            6,
            5,
            &m.clone(),
            6,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn tpmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::tpmv(
            'u',
            'n',
            'n',
            6,
            &m.clone(),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn tpsv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::tpsv(
            'u',
            'n',
            'n',
            6,
            &m.clone(),
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn trmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::trmv(
            'u',
            'n',
            'n',
            6,
            &m.clone(),
            6,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn trsv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    let mut m2: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    for _ in 0..6 {
        m2.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level2::complex::trsv(
            'u',
            'n',
            'n',
            6,
            &m.clone(),
            6,
            &mut m2.clone(),
            1,
        ))
    })
}
