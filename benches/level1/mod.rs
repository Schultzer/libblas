use libblas::level1;
use rand::Rng;
use test::Bencher;

pub mod complex;

#[bench]
pub fn asum(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level1::asum(36, &m.clone(), 1)))
}

#[bench]
fn axpy(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level1::axpy(
            36,
            rng.gen::<f64>(),
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
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level1::copy(36, &m.clone(), -1, &mut m.clone(), -1)))
}

#[bench]
fn dot(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level1::dot(36, &m.clone(), -1, &m.clone(), -1)))
}

#[bench]
fn ddot(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f32> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f32>())
    }
    bh.iter(|| test::black_box(level1::ddot(36, 2.34f32, &m.clone(), -1, &m.clone(), -1)))
}

#[bench]
fn iamax(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| test::black_box(level1::iamax(36, &m.clone(), 1)))
}

#[bench]
fn nrm2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| test::black_box(level1::nrm2(36, &m.clone(), 1)))
}

#[bench]
fn rot(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| {
        test::black_box(level1::rot(
            36,
            &mut m.clone(),
            -1,
            &mut m.clone(),
            -1,
            m[5].clone(),
            m[6].clone(),
        ))
    })
}

#[bench]
fn rotg(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    bh.iter(|| {
        test::black_box(level1::rotg(
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
        ))
    })
}

#[bench]
fn rotm(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| {
        test::black_box(level1::rotm(
            36,
            &mut m.clone(),
            -1,
            &mut m.clone(),
            -1,
            &mut m[0..5],
        ))
    })
}

#[bench]
fn rotmg(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..5 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| {
        test::black_box(level1::rotmg(
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
            &mut rng.gen::<f64>(),
            &mut m.clone(),
        ))
    })
}

#[bench]
fn scal(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }

    bh.iter(|| test::black_box(level1::scal(36, rng.gen::<f64>(), &mut m.clone(), 1)))
}

#[bench]
fn swap(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level1::swap(36, &mut m.clone(), -1, &mut m.clone(), -1)))
}
