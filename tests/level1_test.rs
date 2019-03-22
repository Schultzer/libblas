use libblas::level1;
mod utils;

#[test]
fn iamax() {
    assert_eq!(level1::iamax(6, &vec![1.0, 0.0, 3.0, 4.0, 5.0, 6.0], 1), 6);
    assert_eq!(level1::iamax(3, &vec![1.0, 0.0, 3.0, 4.0, 5.0, 6.0], 2), 3);
    assert_eq!(level1::iamax(0, &vec![1.0, 0.0, 3.0, 4.0, 5.0, 6.0], 0), 0);
    assert_eq!(level1::iamax(1, &vec![1.0, 0.0, 3.0, 4.0, 5.0, 6.0], 1), 1);
}

#[test]
fn asum() {
    assert_eq!(
        level1::asum(6, &vec![1.0, 0.0, 3.0, 4.0, 5.0, 6.0], 1),
        19.0
    );
    assert_eq!(level1::asum(3, &vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2), 9.0);
    assert_eq!(level1::asum(0, &vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2), 0.0);
    assert_eq!(level1::asum(3, &vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 0), 0.0);
}

#[test]
fn axpy() {
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
    let expect = vec![
        28.487472905706035,
        -7.7335447857132165,
        32.14409136136664,
        29.33638278430088,
        9.666040727658336,
        -33.70378597690203,
        -20.896125592422173,
        -8.04363151078942,
        -0.8194978250868761,
        54.861365973632914,
    ];
    level1::axpy(10, 23.0, &x, 1, &mut y, 1);
    assert_eq!(y, expect);

    let x = vec![
        1.2629542848807933098,
        -0.3262333607056494000,
        1.3297992629225006134,
    ];
    let mut y = vec![
        -0.560475646552212603,
        -0.230177489483279957,
        1.558708314149124030,
    ];
    let expect = vec![28.487472905706035, -7.7335447857132165, 32.14409136136664];
    level1::axpy(3, 23.0, &x, 1, &mut y, 1);
    assert_eq!(y, expect);

    let x = vec![
        1.2629542848807933098,
        -0.3262333607056494000,
        1.3297992629225006134,
        1.2724293214294046805,
        0.4146414344564082199,
        -1.5399500419037095433,
        -0.9285670347135380753,
        -0.2947204467905601977,
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
    ];
    let expect = vec![
        28.487472905706035,
        -7.7335447857132165,
        32.14409136136664,
        29.33638278430088,
        9.666040727658336,
        -33.70378597690203,
        -20.896125592422173,
        -8.04363151078942,
    ];
    level1::axpy(8, 23.0, &x, 1, &mut y, 1);
    assert_eq!(y, expect);

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
    let expect = vec![
        28.487472905706035,
        -0.23017748948327996,
        32.14409136136664,
        0.070508391424576,
        9.666040727658336,
        1.715064986883281,
        -20.896125592422173,
        -1.265061234606534,
        -0.8194978250868761,
        -0.44566197009995806,
    ];
    level1::axpy(5, 23.0, &x, 2, &mut y, 2);
    assert_eq!(y, expect);

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
    let expect = vec![
        28.487472905706035,
        -0.23017748948327996,
        32.14409136136664,
        0.070508391424576,
        9.666040727658336,
        1.715064986883281,
        -20.896125592422173,
        -1.265061234606534,
        -0.8194978250868761,
        -0.44566197009995806,
    ];
    level1::axpy(5, 23.0, &x, -2, &mut y, -2);
    assert_eq!(y, expect);
}

#[test]
fn nrm2() {
    assert_eq!(
        level1::nrm2(4, &vec![3.0, 2.0, 3.0, 0.0], 1),
        4.69041575982343
    );
    assert_eq!(level1::nrm2(4, &vec![1.0, 2.0, 3.0, 4.0], 0), 0.0);
    assert_eq!(level1::nrm2(1, &vec![3.0, 2.0, 3.0, 4.0], 1), 3.0);
}

#[test]
fn copy() {
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0];
    level1::copy(4, &x, 1, &mut y, 1);
    assert_eq!(y, vec![1.0, 2.0, 3.0, 4.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    level1::copy(10, &x, 1, &mut y, 1);
    assert_eq!(y, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    level1::copy(7, &x, 1, &mut y, 1);
    assert_eq!(y, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    level1::copy(5, &x, 2, &mut y, 2);
    assert_eq!(y, vec![1.0, 0.0, 3.0, 0.0, 5.0, 0.0, 7.0, 0.0, 9.0, 0.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    level1::copy(5, &x, 2, &mut y, -1);
    assert_eq!(y, vec![9.0, 7.0, 5.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    level1::copy(5, &x, -2, &mut y, 1);
    assert_eq!(y, vec![9.0, 7.0, 5.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0];
    let mut y = vec![0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 9.0, 0.0, 0.0, 0.0];
    level1::copy(0, &x, -2, &mut y, 1);
    assert_eq!(y, vec![0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 9.0, 0.0, 0.0, 0.0]);
}

#[test]
fn dot() {
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let y = vec![5.0, 6.0, 7.0, 8.0];
    assert_eq!(level1::dot(4, &x, 1, &y, 1), 70.0);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 2.0];
    let y = vec![5.0, 6.0, 7.0, 8.0, 5.0, 4.0, 3.0, 3.0, -1.0, -9.0];
    let out =
        (1 * 5 + 2 * 6 + 3 * 7 + 4 * 8 + 5 * 5 + 6 * 4 + 7 * 3 + 8 * 3 + 9 * -1 + 2 * -9) as f64;
    assert_eq!(level1::dot(10, &x, 1, &y, 1), out);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 2.0];
    let y = vec![5.0, 6.0, 7.0, 8.0, 5.0, 4.0, 3.0, 3.0, -1.0, -9.0];
    assert_eq!(level1::dot(0, &x, 1, &y, 1), 0.0);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 2.0];
    let y = vec![5.0, 6.0, 7.0, 8.0, 5.0, 4.0, 3.0, 3.0, -1.0, -9.0];
    let out = (1 * 5 + 3 * 7 + 5 * 5 + 7 * 3 - 9 * 1) as f64;
    assert_eq!(level1::dot(5, &x, 2, &y, 2), out);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 2.0];
    let y = vec![5.0, 6.0, 7.0, 8.0, 5.0, 4.0, 3.0, 3.0, -1.0, -9.0];
    let out = (1 * 5 + 3 * 7 + 5 * 5 + 7 * 3 - 9 * 1) as f64;
    assert_eq!(level1::dot(5, &x, -2, &y, -2), out);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let y = vec![5.0, 6.0, 7.0, 8.0, 5.0, 4.0, 3.0];
    let out = (1 * 5 + 2 * 6 + 3 * 7 + 4 * 8 + 5 * 5 + 6 * 4 + 7 * 3) as f64;
    assert_eq!(level1::dot(7, &x, 1, &y, 1), out);
}

#[test]
fn ddot() {
    let x = vec![1f32, 2f32, 3f32, 4f32];
    let y = vec![5f32, 6f32, 7f32, 8f32];
    assert_eq!(level1::ddot(4, 2f32, &x, 1, &y, 1), 72f64);

    let x = vec![1f32, 2f32, 3f32, 4f32];
    let y = vec![5f32, 6f32, 7f32, 8f32];
    assert_eq!(level1::ddot(0, 2f32, &x, 1, &y, 1), 2f64);

    let x = vec![1f64, 2f64, 3f64, 4f64];
    let y = vec![5f64, 6f64, 7f64, 8f64];
    assert_eq!(level1::ddot(2, 3f64, &x, 2, &y, 2), 29f64);

    let x = vec![1f64, 2f64, 3f64, 4f64];
    let y = vec![5f64, 6f64, 7f64, 8f64];
    assert_eq!(level1::ddot(2, 3f64, &x, -1, &y, -1), 20f64);

    let x = vec![1.0, 2.0, 3.0, 4.0];
    let y = vec![5.0, 6.0, 7.0, 8.0];
    assert_eq!(level1::ddot(2, 3.0, &x, -1, &y, 1), 19f64);
}

#[test]
fn rot() {
    use core::f64::consts::PI;
    let mut x = vec![1.0, 0.0];
    let mut y = vec![0.0, 1.0];
    level1::rot(
        2,
        &mut x,
        1,
        &mut y,
        1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    approximately!(x, vec![0.86602540378443871, 0.5]);
    approximately!(y, vec![-0.5, 0.86602540378443871]);
    // assert_eq!(x, vec![0.86602540378443871, 0.49999999999999994]);
    // assert_eq!(y, vec![-0.49999999999999994, 0.86602540378443871]);

    let mut x = vec![1.0, 0.0];
    let mut y = vec![0.0, 1.0];
    level1::rot(
        2,
        &mut x,
        -1,
        &mut y,
        -1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    approximately!(x, vec![0.86602540378443871, 0.5]);
    approximately!(y, vec![-0.5, 0.86602540378443871]);
    // assert_eq!(x, vec![0.86602540378443871, 0.49999999999999994]);
    // assert_eq!(y, vec![-0.49999999999999994, 0.86602540378443871]);

    let mut x = vec![1.0, 0.0];
    let mut y = vec![0.0, 1.0];
    level1::rot(
        0,
        &mut x,
        -1,
        &mut y,
        -1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    assert_eq!(x, vec![1.0, 0.0]);
    assert_eq!(y, vec![0.0, 1.0]);
}

#[test]
fn rotg() {
    let mut sa = 4.0;
    let mut sb = 2.0;
    let mut c = 0.0;
    let mut s = 0.0;
    level1::rotg(&mut sa, &mut sb, &mut c, &mut s);
    assert_eq!(sa, 4.47213595499958);
    assert_eq!(sb, 0.4472135954999579);
    assert_eq!(c, 0.8944271909999159);
    assert_eq!(s, 0.4472135954999579);

    let mut sa = 0.0;
    let mut sb = 0.0;
    let mut c = 0.0;
    let mut s = 0.0;
    level1::rotg(&mut sa, &mut sb, &mut c, &mut s);
    assert_eq!(sa, 0.0);
    assert_eq!(sb, 0.0);
    assert_eq!(c, 1.0);
    assert_eq!(s, 0.0);

    let mut sa = -2.0;
    let mut sb = -4.0;
    let mut c = 0.0;
    let mut s = 0.0;
    level1::rotg(&mut sa, &mut sb, &mut c, &mut s);
    assert_eq!(sa, -4.47213595499958);
    assert_eq!(sb, 2.2360679774997898);
    assert_eq!(c, 0.4472135954999579);
    assert_eq!(s, 0.8944271909999159);
}

#[test]
fn rotm() {
    use std::f64::NAN;
    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, 1, &mut y, 1, &vec![-2.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(x, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(y, vec![5.0, 6.0, 7.0, 8.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(0, &mut x, 1, &mut y, 1, &vec![-2.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(x, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(y, vec![5.0, 6.0, 7.0, 8.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, 1, &mut y, 1, &vec![-1.0, 2.0, 3.0, 4.0, 5.0]);
    assert_eq!(x, vec![22.0, 28.0, 34.0, 40.0]);
    assert_eq!(y, vec![28.0, 36.0, 44.0, 52.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, 1, &mut y, 1, &vec![0.0, NAN, 3.0, 4.0, NAN]);
    assert_eq!(x, vec![21.0, 26.0, 31.0, 36.0]);
    assert_eq!(y, vec![8.0, 12.0, 16.0, 20.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, 1, &mut y, 1, &vec![1.0, 2.0, NAN, NAN, 3.0]);
    assert_eq!(x, vec![7.0, 10.0, 13.0, 16.0]);
    assert_eq!(y, vec![14.0, 16.0, 18.0, 20.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, -1, &mut y, 1, &vec![1.0, 2.0, NAN, NAN, 3.0]);
    assert_eq!(x, vec![10.0, 11.0, 12.0, 13.0]);
    assert_eq!(y, vec![11.0, 15.0, 19.0, 23.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(4, &mut x, 1, &mut y, -1, &vec![1.0, 2.0, NAN, NAN, 3.0]);
    assert_eq!(x, vec![10.0, 11.0, 12.0, 13.0]);
    assert_eq!(y, vec![11.0, 15.0, 19.0, 23.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(2, &mut x, 2, &mut y, -2, &vec![0.0, NAN, 4.0, 5.0, NAN]);
    assert_eq!(x, vec![36.0, 2.0, 28.0, 4.0]);
    assert_eq!(y, vec![17.0, 6.0, 11.0, 8.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![5.0, 6.0, 7.0, 8.0];
    level1::rotm(2, &mut x, 2, &mut y, -2, &vec![-1.0, 1.0, 2.0, 3.0, 4.0]);
    assert_eq!(x, vec![22.0, 2.0, 18.0, 4.0]);
    assert_eq!(y, vec![26.0, 6.0, 30.0, 8.0]);
}

#[test]
fn rotmg() {
    let mut d1 = -4.0;
    let mut d2 = 2.0;
    let mut x1 = 3.0;
    let mut y1 = 9.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.0);
    assert_eq!(d2, 0.0);
    assert_eq!(x1, 0.0);
    assert_eq!(y1, 9.0);
    assert_eq!(param, [-1.0, 0.0, 0.0, 0.0, 0.0]);

    let mut d1 = 4.0;
    let mut d2 = 0.0;
    let mut x1 = 3.0;
    let mut y1 = 9.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 4.0);
    assert_eq!(d2, 0.0);
    assert_eq!(x1, 3.0);
    assert_eq!(y1, 9.00000000);
    assert_eq!(param, [-2.0, 0.0, 0.0, 0.0, 0.0]);

    let mut d1 = 1.0;
    let mut d2 = 2.0;
    let mut x1 = 3.0;
    let mut y1 = 1.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.81818181818181812);
    assert_eq!(d2, 1.6363636363636362);
    assert_eq!(x1, 3.6666666666666670);
    assert_eq!(y1, 1.0);
    assert_eq!(
        param,
        [0.0, 0.0, -0.33333333333333331, 0.66666666666666663, 0.0]
    );

    let mut d1 = 2.0;
    let mut d2 = -1.0;
    let mut x1 = 3.0;
    let mut y1 = 8.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.0);
    assert_eq!(d2, 0.0);
    assert_eq!(x1, 0.0);
    assert_eq!(y1, 8.0);
    assert_eq!(param, [-1.0, 0.0, 0.0, 0.0, 0.0]);

    let mut d1 = 2.0;
    let mut d2 = 1.0;
    let mut x1 = 3.0;
    let mut y1 = 8.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.78048780487804881);
    assert_eq!(d2, 1.5609756097560976);
    assert_eq!(x1, 10.250000000000000);
    assert_eq!(y1, 8.0000000000000000);
    assert_eq!(param, [1.0, 0.75, 0.0, 0.0, 0.375]);

    let mut d1 = 1f64 / (2 << 23) as f64;
    let mut d2 = 1f64 / (2 << 24) as f64;
    let mut x1 = 3.0;
    let mut y1 = 8.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.39024390243902440);
    assert_eq!(d2, 0.78048780487804881);
    assert_eq!(x1, 2.5024414062500000E-003);
    assert_eq!(y1, 8.0);
    assert_eq!(
        param,
        [
            -1.0,
            1.8310546875000000E-004,
            -2.4414062500000000E-004,
            1.0,
            9.1552734375000000E-005
        ]
    );

    let mut d1 = 1f64 / (2 << 23) as f64;
    let mut d2 = 1f64 / (2 << 24) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.81818181818181812);
    assert_eq!(d2, 0.40909090909090906);
    assert_eq!(x1, 8.9518229166666674E-004);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [
            -1.0,
            2.4414062500000000E-004,
            -2.4414062500000000E-004,
            1.0,
            2.4414062500000000E-004
        ]
    );

    let mut d1 = 2f64 / (2 << 23) as f64;
    let mut d2 = 1f64 / (2 << 23) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 9.7534873268821016E-008);
    assert_eq!(d2, 0.81818181818181812);
    assert_eq!(x1, 3.6666666666666670);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [
            -1.0,
            1.0,
            -1.6276041666666666E-004,
            0.33333333333333331,
            2.4414062500000000E-004
        ]
    );

    let mut d1 = (2 << 23) as f64;
    let mut d2 = 2f64 * (2 << 23) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 8882055.5294117648);
    assert_eq!(d2, 1.0588235294117647);
    assert_eq!(x1, 5.6666666666666661);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [-1.0, 1.0, -2730.6666666666665, 1.3333333333333333, 4096.0]
    );

    let mut d1 = 2f64 * (2 << 23) as f64;
    let mut d2 = 1f64 * (2 << 23) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    level1::rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 1.6363636363636362);
    assert_eq!(d2, 13726813.090909090);
    assert_eq!(x1, 15018.666666666668);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [-1.0, 4096.0, -0.66666666666666663, 1365.3333333333333, 1.0]
    );
}

#[test]
fn scal() {
    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    level1::scal(4, 1.0, &mut x, 1);
    assert_eq!(x, vec![1.0, 2.0, 3.0, 4.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    level1::scal(4, 2.0, &mut x, 1);
    assert_eq!(x, vec![2.0, 4.0, 6.0, 8.0]);

    let mut x = vec![3.0, 1.0, 2.0, 3.0, 4.0];
    level1::scal(5, 2.0, &mut x, 1);
    assert_eq!(x, vec![6.0, 2.0, 4.0, 6.0, 8.0]);

    let mut x = vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0];
    level1::scal(7, 2.0, &mut x, 1);
    assert_eq!(x, vec![-2.0, 4.0, 6.0, 2.0, 4.0, 6.0, 8.0]);

    let mut x = vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0];
    level1::scal(3, 2.0, &mut x, 2);
    assert_eq!(x, vec![-2.0, 2.0, 6.0, 1.0, 4.0, 3.0, 4.0]);

    let mut x = vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0];
    level1::scal(6, 2.0, &mut x, 0);
    assert_eq!(x, vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0]);

    let mut x = vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0];
    level1::scal(0, 2.0, &mut x, 2);
    assert_eq!(x, vec![-1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn swap() {
    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![1.0, 2.0, 3.0, 4.0];
    level1::swap(4, &mut x, 1, &mut y, -1);
    assert_eq!(x, vec![4.0, 3.0, 2.0, 1.0]);
    assert_eq!(y, vec![4.0, 3.0, 2.0, 1.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![1.0, 2.0, 3.0, 4.0];
    level1::swap(4, &mut x, -1, &mut y, 1);
    assert_eq!(x, vec![4.0, 3.0, 2.0, 1.0]);
    assert_eq!(y, vec![4.0, 3.0, 2.0, 1.0]);

    let mut x = vec![1.0, 2.0, 3.0, 4.0];
    let mut y = vec![7.0, 8.0, 9.0, 5.0];
    level1::swap(4, &mut x, 1, &mut y, 1);
    assert_eq!(x, vec![7.0, 8.0, 9.0, 5.0]);
    assert_eq!(y, vec![1.0, 2.0, 3.0, 4.0]);

    let mut x = vec![1.0, 2.0, 3.0];
    let mut y = vec![7.0, 8.0, 9.0];
    level1::swap(3, &mut x, 1, &mut y, 1);
    assert_eq!(x, vec![7.0, 8.0, 9.0]);
    assert_eq!(y, vec![1.0, 2.0, 3.0]);

    let mut x = vec![1.0, 2.0];
    let mut y = vec![7.0, 8.0];
    level1::swap(2, &mut x, 1, &mut y, 1);
    assert_eq!(x, vec![7.0, 8.0]);
    assert_eq!(y, vec![1.0, 2.0]);

    let mut x = vec![1.0, 2.0];
    let mut y = vec![7.0, 8.0];
    level1::swap(0, &mut x, 1, &mut y, 1);
    assert_eq!(x, vec![1.0, 2.0]);
    assert_eq!(y, vec![7.0, 8.0]);
}
