use libblas::level1;
use num_complex::Complex;
use rand::Rng;
use test::Bencher;

#[bench]
pub fn asum(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| test::black_box(level1::complex::asum(36, &m.clone(), 1)))
}

#[bench]
fn axpy(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level1::complex::axpy(
            36,
            &Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            &m.clone(),
            -1,
            &mut m.clone(),
            -1,
        ))
    })
}

#[bench]
fn copy(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level1::complex::copy(
            36,
            &m.clone(),
            -1,
            &mut m.clone(),
            -1,
        ))
    })
}

#[bench]
fn dotc(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| test::black_box(level1::complex::dotc(36, &m.clone(), -1, &m.clone(), -1)))
}

#[bench]
fn dotu(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| test::black_box(level1::complex::dotu(36, &m.clone(), -1, &m.clone(), -1)))
}

#[bench]
fn iamax(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| test::black_box(level1::complex::iamax(36, &m.clone(), 1)))
}

#[bench]
fn nrm2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| test::black_box(level1::complex::nrm2(36, &m.clone(), 1)))
}

#[bench]
fn rot(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level1::complex::rot(
            36,
            &mut m.clone(),
            -1,
            &mut m.clone(),
            -1,
            m[5].re.clone(),
            m[6].re.clone(),
        ))
    })
}

#[bench]
fn rotg(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    bh.iter(|| {
        test::black_box(level1::complex::rotg(
            &mut Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            &mut Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            &mut rng.gen::<f64>(),
            &mut Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
        ))
    })
}

#[bench]
fn scal(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }

    bh.iter(|| {
        test::black_box(level1::complex::scal(
            36,
            Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            &mut m.clone(),
            1,
        ))
    })
}

#[bench]
fn sscal(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }

    bh.iter(|| {
        test::black_box(level1::complex::sscal(
            36,
            rng.gen::<f64>(),
            &mut m.clone(),
            1,
        ))
    })
}

#[bench]
fn swap(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<Complex<f64>> = Vec::new();
    for _ in 0..36 {
        m.push(Complex::new(rng.gen::<f64>(), rng.gen::<f64>()))
    }
    bh.iter(|| {
        test::black_box(level1::complex::swap(
            36,
            &mut m.clone(),
            -1,
            &mut m.clone(),
            -1,
        ))
    })
}
