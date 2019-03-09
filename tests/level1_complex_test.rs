use blasrs::level1::complex;
use num_complex::Complex;
mod utils;

#[test]
fn axpy() {
    let x = vec![
        Complex::new(1.2629542848807933098, 0.76359346114045956),
        Complex::new(-0.3262333607056494000, -0.79900924898936820),
        Complex::new(1.3297992629225006134, -1.14765700923635139),
        Complex::new(1.2724293214294046805, -0.28946157368822334),
        Complex::new(0.4146414344564082199, -0.29921511789731614),
        Complex::new(-1.5399500419037095433, -0.41151083279506701),
        Complex::new(-0.9285670347135380753, 0.25222344815613229),
        Complex::new(-0.2947204467905601977, -0.89192112728456863),
        Complex::new(-0.0057671727475369552, 0.43568329935571865),
        Complex::new(2.4046533888579508798, -1.23753842192995811),
    ];
    let mut y = vec![
        Complex::new(-0.560475646552212603, 1.22408179743946155),
        Complex::new(-0.230177489483279957, 0.35981382705736381),
        Complex::new(1.558708314149124030, 0.40077145059405217),
        Complex::new(0.070508391424576003, 0.11068271594511971),
        Complex::new(0.129287735160946243, -0.55584113475407493),
        Complex::new(1.715064986883281017, 1.78691313680307817),
        Complex::new(0.460916205989202299, 0.49785047822923939),
        Complex::new(-1.265061234606533969, -1.96661715662963821),
        Complex::new(-0.686852851893526073, 0.70135590156368555),
        Complex::new(-0.445661970099958060, -0.47279140772793404),
    ];
    let a = Complex::new(23.0, 32.0);
    let expect = vec![
        Complex::new(4.052482149211329, 59.2012685198554109),
        Complex::new(17.834751181946569, -28.4568664422788906),
        Complex::new(68.869115656929893, 16.5582366516779906),
        Complex::new(38.599153142324027, 34.1708048068569283),
        Complex::new(19.240924500372454, 5.8307370562127172),
        Complex::new(-20.535439327459894, -56.9562373584021628),
        Complex::new(-28.967275933418406, -23.4151553250129361),
        Complex::new(20.497844562316779, -31.9118573814726432),
        Complex::new(-14.761363404469874, 10.5375222588240316),
        Complex::new(94.462595475391566, 48.0127333313374578),
    ];

    complex::axpy(10, &a, &x, 1, &mut y, 1);
    assert_eq!(y, expect);

    // let x = vec![
    //         Complex::new(1.2629542848807933098,0.76359346114045956),
    //         Complex::new(-0.3262333607056494000, -0.79900924898936820),
    //         Complex::new(1.3297992629225006134, -1.14765700923635139),
    //         Complex::new(1.2724293214294046805, -0.28946157368822334),
    //         Complex::new(0.4146414344564082199, -0.29921511789731614),
    //         Complex::new(-1.5399500419037095433, -0.41151083279506701),
    //         Complex::new(-0.9285670347135380753,  0.25222344815613229),
    //         Complex::new(-0.2947204467905601977, -0.89192112728456863),
    //         Complex::new(-0.0057671727475369552, 0.43568329935571865),
    //         Complex::new(2.4046533888579508798, -1.23753842192995811)
    //     ],
    // };
    let mut y = vec![
        Complex::new(-0.560475646552212603, 1.22408179743946155),
        Complex::new(-0.230177489483279957, 0.35981382705736381),
        Complex::new(1.558708314149124030, 0.40077145059405217),
        Complex::new(0.070508391424576003, 0.11068271594511971),
        Complex::new(0.129287735160946243, -0.55584113475407493),
        Complex::new(1.715064986883281017, 1.78691313680307817),
        Complex::new(0.460916205989202299, 0.49785047822923939),
        Complex::new(-1.265061234606533969, -1.96661715662963821),
        Complex::new(-0.686852851893526073, 0.70135590156368555),
        Complex::new(-0.445661970099958060, -0.47279140772793404),
    ];
    let expect = vec![
        Complex::new(4.052482149211329, 59.2012685198554109),
        Complex::new(-0.230177489483279957, 0.35981382705736381),
        Complex::new(68.869115656929893, 16.5582366516779906),
        Complex::new(0.070508391424576003, 0.11068271594511971),
        Complex::new(19.240924500372454, 5.8307370562127172),
        Complex::new(1.715064986883281017, 1.78691313680307817),
        Complex::new(-28.967275933418406, -23.4151553250129361),
        Complex::new(-1.265061234606533969, -1.96661715662963821),
        Complex::new(-14.761363404469874, 10.5375222588240316),
        Complex::new(-0.445661970099958060, -0.47279140772793404),
    ];

    complex::axpy(5, &a, &x, 2, &mut y, 2);
    assert_eq!(y, expect);

    let x = vec![Complex::new(3.0, 4.0)];
    let mut y = vec![Complex::new(1.0, 2.0)];
    let a = Complex::new(0.0, 0.0);
    complex::axpy(1, &a, &x, 1, &mut y, 1);
    assert_eq!(y, vec![Complex::new(1.0, 2.0)]);

    let x = vec![Complex::new(9.0, 4.0)];
    let mut y = vec![Complex::new(3.0, 4.0)];
    let a = Complex::new(1.0, 2.0);
    complex::axpy(0, &a, &x, 1, &mut y, 1);
    assert_eq!(y, vec![Complex::new(3.0, 4.0)]);

    let x = vec![
        Complex::new(1.2629542848807933098, 0.76359346114045956),
        Complex::new(-0.3262333607056494000, -0.79900924898936820),
        Complex::new(1.3297992629225006134, -1.14765700923635139),
        Complex::new(1.2724293214294046805, -0.28946157368822334),
        Complex::new(0.4146414344564082199, -0.29921511789731614),
        Complex::new(-1.5399500419037095433, -0.41151083279506701),
        Complex::new(-0.9285670347135380753, 0.25222344815613229),
        Complex::new(-0.2947204467905601977, -0.89192112728456863),
        Complex::new(-0.0057671727475369552, 0.43568329935571865),
        Complex::new(2.4046533888579508798, -1.23753842192995811),
    ];
    let mut y = vec![
        Complex::new(-0.560475646552212603, 1.22408179743946155),
        Complex::new(-0.230177489483279957, 0.35981382705736381),
        Complex::new(1.558708314149124030, 0.40077145059405217),
        Complex::new(0.070508391424576003, 0.11068271594511971),
        Complex::new(0.129287735160946243, -0.55584113475407493),
        Complex::new(1.715064986883281017, 1.78691313680307817),
        Complex::new(0.460916205989202299, 0.49785047822923939),
        Complex::new(-1.265061234606533969, -1.96661715662963821),
        Complex::new(-0.686852851893526073, 0.70135590156368555),
        Complex::new(-0.445661970099958060, -0.47279140772793404),
    ];

    let a = Complex::new(23.0, 32.0);

    let expect = vec![
        Complex::new(4.052482149211329, 59.2012685198554109),
        Complex::new(-0.230177489483279957, 0.35981382705736381),
        Complex::new(68.869115656929893, 16.5582366516779906),
        Complex::new(0.070508391424576003, 0.11068271594511971),
        Complex::new(19.240924500372454, 5.8307370562127172),
        Complex::new(1.715064986883281017, 1.78691313680307817),
        Complex::new(-28.967275933418406, -23.4151553250129361),
        Complex::new(-1.265061234606533969, -1.96661715662963821),
        Complex::new(-14.761363404469874, 10.5375222588240316),
        Complex::new(-0.445661970099958060, -0.47279140772793404),
    ];

    complex::axpy(5, &a, &x, -2, &mut y, -2);
    assert_eq!(y, expect);
}

#[test]
fn copy() {
    let x = vec![
        Complex::new(1.0, 5.0),
        Complex::new(2.0, 6.0),
        Complex::new(3.0, 7.0),
        Complex::new(4.0, 8.0),
    ];
    let mut y = vec![
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];
    let expect = vec![
        Complex::new(1.0, 5.0),
        Complex::new(2.0, 6.0),
        Complex::new(3.0, 7.0),
        Complex::new(4.0, 8.0),
    ];
    complex::copy(4, &x, 1, &mut y, 1);
    assert_eq!(y, expect);

    let mut y = vec![
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];
    complex::copy(4, &x, -1, &mut y, -1);
    assert_eq!(y, expect);

    let mut y = vec![
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];
    let expect = vec![
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];
    complex::copy(0, &x, -1, &mut y, -1);
    assert_eq!(y, expect);
}

#[test]
fn dotc() {
    let x = vec![
        Complex::new(1.0, 5.0),
        Complex::new(2.0, 6.0),
        Complex::new(3.0, 7.0),
        Complex::new(4.0, 8.0),
    ];
    let y = vec![
        Complex::new(5.0, 9.0),
        Complex::new(6.0, 10.0),
        Complex::new(7.0, 11.0),
        Complex::new(8.0, 12.0),
    ];
    assert_eq!(complex::dotc(4, &x, 1, &y, 1), Complex::new(348.0, -64.0));
    assert_eq!(complex::dotc(4, &x, -1, &y, -1), Complex::new(348.0, -64.0));

    let x = vec![Complex::new(1.0, 5.0)];
    let y = vec![Complex::new(5.0, 9.0)];
    assert_eq!(complex::dotc(0, &x, -1, &y, -1), Complex::new(0.0, 0.0));
}

#[test]
fn dotu() {
    let x = vec![
        Complex::new(1.0, 5.0),
        Complex::new(2.0, 6.0),
        Complex::new(3.0, 7.0),
        Complex::new(4.0, 8.0),
    ];
    let y = vec![
        Complex::new(5.0, 9.0),
        Complex::new(6.0, 10.0),
        Complex::new(7.0, 11.0),
        Complex::new(8.0, 12.0),
    ];
    assert_eq!(complex::dotu(4, &x, 1, &y, 1), Complex::new(-208.0, 284.0));
    assert_eq!(
        complex::dotu(4, &x, -1, &y, -1),
        Complex::new(-208.0, 284.0)
    );

    let x = vec![Complex::new(1.0, 5.0)];
    let y = vec![Complex::new(5.0, 9.0)];
    assert_eq!(complex::dotu(0, &x, -1, &y, -1), Complex::new(0.0, 0.0));
}

#[test]
fn rot() {
    use core::f64::consts::PI;
    let mut x = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
    let mut y = vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)];
    complex::rot(
        2,
        &mut x,
        1,
        &mut y,
        1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    assert_eq!(
        x,
        vec![
            Complex::new(0.8660254037844387, 0.0),
            Complex::new(0.49999999999999994, 0.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(-0.49999999999999994, 0.0),
            Complex::new(0.8660254037844387, 0.0)
        ]
    );
    // assert_eq!(x, vec![Complex::new(0.8660254037844387, 0.0), Complex::new(0.500000000000000, 0.0)]);
    // assert_eq!(y, vec![Complex::new(-0.500000000000000, 0.0), Complex::new(0.8660254037844387, 0.0)]);

    let mut x = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
    let mut y = vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)];
    complex::rot(
        2,
        &mut x,
        -1,
        &mut y,
        -1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    assert_eq!(
        x,
        vec![
            Complex::new(0.8660254037844387, 0.0),
            Complex::new(0.49999999999999994, 0.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(-0.49999999999999994, 0.0),
            Complex::new(0.8660254037844387, 0.0)
        ]
    );
    // assert_eq!(x, vec![Complex::new(0.8660254037844387, 0.0), Complex::new(0.500000000000000, 0.0)]);
    // assert_eq!(y, vec![Complex::new(-0.500000000000000, 0.0), Complex::new(0.8660254037844387, 0.0)]);

    let mut x = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
    let mut y = vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)];
    complex::rot(
        0,
        &mut x,
        -1,
        &mut y,
        -1,
        (PI * (1.0 / 6.0)).cos(),
        (PI * (1.0 / 6.0)).sin(),
    );
    assert_eq!(x, vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);
    assert_eq!(y, vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]);
}

#[test]
fn rotg() {
    let mut a = Complex::new(11.0, 19.0);
    let mut b = Complex::new(34.0, 23.0);
    let mut c = 0.0;
    let mut s = Complex::new(0.0, 0.0);
    complex::rotg(&mut a, &mut b, &mut c, &mut s);
    assert_eq!(a, Complex::new(23.323763103564644, 40.286499906157111));
    assert_eq!(b, Complex::new(34.000000, 23.0000000));
    assert_eq!(c, 0.47162200847078728);
    assert_eq!(s, Complex::new(0.79353827566350310, 0.38453827661622281));

    let mut a = Complex::new(0.0, 0.0);
    let mut b = Complex::new(-1.0, -1.0);
    let mut c = 0.0;
    let mut s = Complex::new(0.0, 0.0);
    complex::rotg(&mut a, &mut b, &mut c, &mut s);
    assert_eq!(a, Complex::new(-1.0, -1.0));
    assert_eq!(b, Complex::new(-1.0, -1.0));
    assert_eq!(c, 0.0);
    assert_eq!(s, Complex::new(1.0, 0.0));
}

#[test]
fn sscal() {
    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    complex::sscal(6, 2.0, &mut x, 1);
    assert_eq!(
        x,
        vec![
            Complex::new(2.0, 14.0),
            Complex::new(4.0, 16.0),
            Complex::new(6.0, 18.0),
            Complex::new(8.0, 20.0),
            Complex::new(10.0, 22.0),
            Complex::new(12.0, 24.0)
        ]
    );

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    complex::sscal(3, 2.0, &mut x, 2);
    assert_eq!(
        x,
        vec![
            Complex::new(2.0, 14.0),
            Complex::new(2.0, 8.0),
            Complex::new(6.0, 18.0),
            Complex::new(4.0, 10.0),
            Complex::new(10.0, 22.0),
            Complex::new(6.0, 12.0)
        ]
    );

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    let expect = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    complex::sscal(3, 2.0, &mut x, 0);
    assert_eq!(x, expect);

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    complex::sscal(0, 2.0, &mut x, 1);
    assert_eq!(x, expect);

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    complex::sscal(6, 0.0, &mut x, 1);
    assert_eq!(
        x,
        vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0)
        ]
    );
}

#[test]
fn scal() {
    let mut x = vec![
        Complex::new(1.26295428488079331, -0.9285670347135380753),
        Complex::new(-0.32623336070564940, -0.2947204467905601977),
        Complex::new(1.32979926292250061, -0.0057671727475369552),
        Complex::new(1.27242932142940468, 2.4046533888579508798),
        Complex::new(0.41464143445640822, 0.7635934611404595618),
        Complex::new(-1.53995004190370954, -0.7990092489893682037),
    ];
    let a = Complex::new(2.0, 4.0);
    complex::scal(6, a, &mut x, 1);
    assert_eq!(
        x,
        vec![
            Complex::new(6.24017670861573848, 3.1946830700960973),
            Complex::new(0.52641506575094199, -1.8943743364037180),
            Complex::new(2.68266721683514886, 5.3076627061949289),
            Complex::new(-7.07375491257299416, 9.8990240634335205),
            Complex::new(-2.22509097564902181, 3.1857526601065520),
            Complex::new(0.11613691215005373, -7.7578186655935744)
        ]
    );

    let mut x = vec![
        Complex::new(1.26295428488079331, -0.9285670347135380753),
        Complex::new(-0.32623336070564940, -0.2947204467905601977),
        Complex::new(1.32979926292250061, -0.0057671727475369552),
        Complex::new(1.27242932142940468, 2.4046533888579508798),
        Complex::new(0.41464143445640822, 0.7635934611404595618),
        Complex::new(-1.53995004190370954, -0.7990092489893682037),
    ];
    complex::scal(3, a, &mut x, 2);
    assert_eq!(
        x,
        vec![
            Complex::new(6.24017670861573848, 3.1946830700960973),
            Complex::new(-0.32623336070564940, -0.2947204467905601977),
            Complex::new(2.68266721683514886, 5.3076627061949289),
            Complex::new(1.27242932142940468, 2.4046533888579508798),
            Complex::new(-2.22509097564902181, 3.1857526601065520),
            Complex::new(-1.53995004190370954, -0.7990092489893682037)
        ]
    );

    let mut x = vec![
        Complex::new(1.26295428488079331, -0.9285670347135380753),
        Complex::new(-0.32623336070564940, -0.2947204467905601977),
        Complex::new(1.32979926292250061, -0.0057671727475369552),
        Complex::new(1.27242932142940468, 2.4046533888579508798),
        Complex::new(0.41464143445640822, 0.7635934611404595618),
        Complex::new(-1.53995004190370954, -0.7990092489893682037),
    ];
    complex::scal(0, a, &mut x, 0);
    assert_eq!(
        x,
        vec![
            Complex::new(1.26295428488079331, -0.9285670347135380753),
            Complex::new(-0.32623336070564940, -0.2947204467905601977),
            Complex::new(1.32979926292250061, -0.0057671727475369552),
            Complex::new(1.27242932142940468, 2.4046533888579508798),
            Complex::new(0.41464143445640822, 0.7635934611404595618),
            Complex::new(-1.53995004190370954, -0.7990092489893682037)
        ]
    );

    let mut x = vec![
        Complex::new(1.26295428488079331, -0.9285670347135380753),
        Complex::new(-0.32623336070564940, -0.2947204467905601977),
        Complex::new(1.32979926292250061, -0.0057671727475369552),
        Complex::new(1.27242932142940468, 2.4046533888579508798),
        Complex::new(0.41464143445640822, 0.7635934611404595618),
        Complex::new(-1.53995004190370954, -0.7990092489893682037),
    ];
    complex::scal(3, a, &mut x, 0);
    assert_eq!(
        x,
        vec![
            Complex::new(1.26295428488079331, -0.9285670347135380753),
            Complex::new(-0.32623336070564940, -0.2947204467905601977),
            Complex::new(1.32979926292250061, -0.0057671727475369552),
            Complex::new(1.27242932142940468, 2.4046533888579508798),
            Complex::new(0.41464143445640822, 0.7635934611404595618),
            Complex::new(-1.53995004190370954, -0.7990092489893682037)
        ]
    );
}

#[test]
fn swap() {
    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    let mut y = vec![
        Complex::new(13.0, 74.0),
        Complex::new(23.0, 84.0),
        Complex::new(33.0, 94.0),
        Complex::new(43.0, 104.0),
        Complex::new(53.0, 114.0),
        Complex::new(63.0, 124.0),
    ];

    complex::swap(6, &mut x, 1, &mut y, 1);
    assert_eq!(
        x,
        vec![
            Complex::new(13.0, 74.0),
            Complex::new(23.0, 84.0),
            Complex::new(33.0, 94.0),
            Complex::new(43.0, 104.0),
            Complex::new(53.0, 114.0),
            Complex::new(63.0, 124.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(1.0, 7.0),
            Complex::new(2.0, 8.0),
            Complex::new(3.0, 9.0),
            Complex::new(4.0, 10.0),
            Complex::new(5.0, 11.0),
            Complex::new(6.0, 12.0)
        ]
    );

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    let mut y = vec![
        Complex::new(13.0, 74.0),
        Complex::new(23.0, 84.0),
        Complex::new(33.0, 94.0),
        Complex::new(43.0, 104.0),
        Complex::new(53.0, 114.0),
        Complex::new(63.0, 124.0),
    ];

    complex::swap(6, &mut x, -1, &mut y, 1);
    assert_eq!(
        x,
        vec![
            Complex::new(63.0, 124.0),
            Complex::new(53.0, 114.0),
            Complex::new(43.0, 104.0),
            Complex::new(33.0, 94.0),
            Complex::new(23.0, 84.0),
            Complex::new(13.0, 74.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(6.0, 12.0),
            Complex::new(5.0, 11.0),
            Complex::new(4.0, 10.0),
            Complex::new(3.0, 9.0),
            Complex::new(2.0, 8.0),
            Complex::new(1.0, 7.0)
        ]
    );

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    let mut y = vec![
        Complex::new(13.0, 74.0),
        Complex::new(23.0, 84.0),
        Complex::new(33.0, 94.0),
        Complex::new(43.0, 104.0),
        Complex::new(53.0, 114.0),
        Complex::new(63.0, 124.0),
    ];

    complex::swap(6, &mut x, 1, &mut y, -1);
    assert_eq!(
        x,
        vec![
            Complex::new(63.0, 124.0),
            Complex::new(53.0, 114.0),
            Complex::new(43.0, 104.0),
            Complex::new(33.0, 94.0),
            Complex::new(23.0, 84.0),
            Complex::new(13.0, 74.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(6.0, 12.0),
            Complex::new(5.0, 11.0),
            Complex::new(4.0, 10.0),
            Complex::new(3.0, 9.0),
            Complex::new(2.0, 8.0),
            Complex::new(1.0, 7.0)
        ]
    );

    let mut x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    let mut y = vec![
        Complex::new(13.0, 74.0),
        Complex::new(23.0, 84.0),
        Complex::new(33.0, 94.0),
        Complex::new(43.0, 104.0),
        Complex::new(53.0, 114.0),
        Complex::new(63.0, 124.0),
    ];
    complex::swap(0, &mut x, 1, &mut y, -1);
    assert_eq!(
        x,
        vec![
            Complex::new(1.0, 7.0),
            Complex::new(2.0, 8.0),
            Complex::new(3.0, 9.0),
            Complex::new(4.0, 10.0),
            Complex::new(5.0, 11.0),
            Complex::new(6.0, 12.0)
        ]
    );
    assert_eq!(
        y,
        vec![
            Complex::new(13.0, 74.0),
            Complex::new(23.0, 84.0),
            Complex::new(33.0, 94.0),
            Complex::new(43.0, 104.0),
            Complex::new(53.0, 114.0),
            Complex::new(63.0, 124.0)
        ]
    );
}

#[test]
fn iamax() {
    let x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(0.0, 0.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    assert_eq!(complex::iamax(6, &x, 1), 6);

    let x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(2.0, 8.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    assert_eq!(complex::iamax(3, &x, 2), 3);
    assert_eq!(complex::iamax(6, &x, 0), 0);
    assert_eq!(complex::iamax(0, &x, 1), 0);

    let x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(0.0, 1.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    assert_eq!(complex::iamax(1, &x, 1), 1);
}

#[test]
fn asum() {
    let x = vec![
        Complex::new(1.0, 7.0),
        Complex::new(0.0, 0.0),
        Complex::new(3.0, 9.0),
        Complex::new(4.0, 10.0),
        Complex::new(5.0, 11.0),
        Complex::new(6.0, 12.0),
    ];
    assert_eq!(complex::asum(6, &x, 1), 68.0);
    assert_eq!(complex::asum(3, &x, 2), 36.0);
    assert_eq!(complex::asum(0, &x, 2), 0.0);
    assert_eq!(complex::asum(0, &x, 0), 0.0);
}

#[test]
fn nrm2() {
    let x = vec![
        Complex::new(0.0, 0.76359346114045956),
        Complex::new(0.0, 0.0),
        Complex::new(1.3297992629225006134, -1.14765700923635139),
        Complex::new(1.2724293214294046805, -0.28946157368822334),
        Complex::new(0.4146414344564082199, -0.29921511789731614),
        Complex::new(-1.5399500419037095433, -0.41151083279506701),
        Complex::new(-0.9285670347135380753, 0.25222344815613229),
        Complex::new(-0.2947204467905601977, -0.89192112728456863),
        Complex::new(-0.0057671727475369552, 0.43568329935571865),
        Complex::new(2.4046533888579508798, -1.23753842192995811),
    ];
    assert_approx!(complex::nrm2(10, &x, 1), 4.1815805452999522);
    assert_eq!(complex::nrm2(0, &x, 1), 0.0);
    assert_eq!(complex::nrm2(10, &x, 0), 0.0);
}
