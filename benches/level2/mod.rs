use libblas::level2;
use rand::Rng;
use test::Bencher;

pub mod complex;

#[bench]
fn gbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::gbmv(
            'n',
            6,
            6,
            1,
            1,
            1.5,
            &m.clone(),
            6,
            &m2.clone(),
            1,
            2.5,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn gemv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::gemv(
            'n',
            6,
            6,
            1.5,
            &m.clone(),
            6,
            &m2.clone(),
            1,
            2.5,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn ger(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::ger(
            6,
            6,
            1.5,
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
fn sbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::sbmv(
            'u',
            6,
            5,
            1.5,
            &m.clone(),
            6,
            &m2.clone(),
            1,
            2.5,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn spmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..64 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..8 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::spmv(
            'u',
            8,
            1.0,
            &m,
            &m2.clone(),
            1,
            0.25,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn spr(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..64 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..8 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level2::spr('u', 8, 1.0, &m2.clone(), 1, &mut m.clone())))
}

#[bench]
fn spr2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..64 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..8 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::spr2(
            'u',
            8,
            1.0,
            &m2.clone(),
            1,
            &m2.clone(),
            1,
            &mut m.clone(),
        ))
    })
}

#[bench]
fn symv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::symv(
            'u',
            6,
            1.0,
            &m.clone(),
            6,
            &m2.clone(),
            1,
            1.0,
            &mut m2.clone(),
            1,
        ))
    })
}

#[bench]
fn syr(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| test::black_box(level2::syr('u', 6, 1.0, &m2.clone(), 1, &mut m.clone(), 6)))
}

#[bench]
fn syr2(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::syr2(
            'u',
            6,
            1.0,
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
fn tbmv(bh: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::tbmv(
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
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::tbsv(
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
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::tpmv(
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
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::tpsv(
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
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::trmv(
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
    let mut m: Vec<f64> = Vec::new();
    let mut m2: Vec<f64> = Vec::new();
    for _ in 0..36 {
        m.push(rng.gen::<f64>())
    }
    for _ in 0..6 {
        m2.push(rng.gen::<f64>())
    }
    bh.iter(|| {
        test::black_box(level2::trsv(
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
