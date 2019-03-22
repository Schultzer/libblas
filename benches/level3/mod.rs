use libblas::level3;
use rand::Rng;
use test::Bencher;

pub mod complex;

#[bench]
fn gemm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::gemm(
            'n',
            'n',
            6,
            6,
            4,
            1.5,
            &m.clone(),
            6,
            &m.clone(),
            6,
            2.5,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn symm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::symm(
            'l',
            'u',
            6,
            6,
            1.5,
            &m.clone(),
            6,
            &m.clone(),
            6,
            2.5,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn syrk(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::syrk(
            'l',
            't',
            6,
            3,
            0.2,
            &m.clone(),
            6,
            1.0,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn syr2k(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::syr2k(
            'l',
            't',
            6,
            3,
            0.2,
            &m.clone(),
            6,
            &m.clone(),
            6,
            1.0,
            &mut m.clone(),
            6,
        ))
    })
}

#[bench]
fn trmm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::trmm(
            'l',
            'u',
            'n',
            'u',
            3,
            0,
            0.25,
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
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level3::trsm(
            'l',
            'u',
            'n',
            'u',
            3,
            0,
            0.25,
            &m.clone(),
            6,
            &mut m.clone(),
            6,
        ))
    })
}
