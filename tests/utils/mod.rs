#[macro_export]
macro_rules! assert_approx {
    ( $( $l:expr, $r:expr ),* ) => {
        {
            $(
                let (a, b) = (&$l, &$r);
                assert!(((*a - *b) as f64).abs() <= 1.0e-6, "{} is not approximately equal to {}", *a, *b);
            )*
        }
    };
    ( $( $l:expr, $r:expr, $w:expr),* ) => {
        {
            $(
                let (a, b) = (&$l, &$r);
                assert!(((*a - *b) as f64).abs() <= $w, "{} is not approximately equal to {}", *a, *b);
            )*
        }
    };
}

#[macro_export]
macro_rules! approximately {
    ( $( $l:expr, $r:expr ),* ) => {
        {
            $(
                for i in 0..$l.len() {
                    let (a, b) = (&$l[i], &$r[i]);
                    assert!(((*a - *b) as f64).abs() <= 1.0e-6, "{} is not approximately equal to {}", *a, *b);
                }
            )*
        }
    };
    ( $( $l:expr, $r:expr, $w:expr),* ) => {
        {
            $(
                for i in 0..$l.len() {
                    let (a, b) = (&$l[i], &$r[i]);
                        assert!(((*a - *b) as f64).abs() <= $w, "{} is not approximately equal to {}", *a, *b);
                }
            )*
        }
    };
}

#[macro_export]
macro_rules! capproximately {
    ( $( $l:expr, $r:expr ),* ) => {
        {
            $(
                for i in 0..$l.len() {
                    let (a, b) = (&$l[i], &$r[i]);
                    assert!(((a.re - b.re) as f64).abs() <= 1.0e-6, "{} is not approximately equal to {}", a.re, b.re);
                    assert!(((a.im - b.im) as f64).abs() <= 1.0e-6, "{} is not approximately equal to {}", a.im, b.im);
                }
            )*
        }
    };
    ( $( $l:expr, $r:expr, $w:expr),* ) => {
        {
            $(
                for i in 0..$l.len() {
                    let (a, b) = (&$l[i], &$r[i]);
                    println!("{:?} {:?}", a, b);
                    assert!(((a.re - b.re) as f64).abs() <= $w, "{} is not approximately equal to {}", a.re, b.re);
                    assert!(((a.im - b.im) as f64).abs() <= $w, "{} is not approximately equal to {}", a.im, b.im);
                }
            )*
        }
    };
}

#[macro_export]
macro_rules! assert_approx_eq {
    ($a:expr, $b:expr) => {{
        let eps = 1.0e-6;
        let (a, b) = (&$a, &$b);
        assert!(
            ((*a - *b) as f64).abs() < eps,
            "assertion failed: `(left !== right)` \
             (left: `{:?}`, right: `{:?}`, expect diff: `{:?}`, real diff: `{:?}`)",
            *a,
            *b,
            eps,
            ((*a - *b) as f64).abs()
        );
    }};
    ($a:expr, $b:expr, $eps:expr) => {{
        let (a, b) = (&$a, &$b);
        let eps = $eps;
        assert!(
            ((*a - *b) as f64).abs() < eps,
            "assertion failed: `(left !== right)` \
             (left: `{:?}`, right: `{:?}`, expect diff: `{:?}`, real diff: `{:?}`)",
            *a,
            *b,
            eps,
            ((*a - *b) as f64).abs()
        );
    }};
}

#[macro_export]
macro_rules! assert_approx_eq_cplx {
    ($a:expr, $b:expr) => {{
        let eps = 1.0e-6;
        let (a, b) = (&$a, &$b);
        assert!(
            ((a.re - b.re) as f64).abs() <= eps,
            "{} is not approximately equal to {}",
            a.re,
            b.re
        );
        assert!(
            ((a.im - b.im) as f64).abs() <= eps,
            "{} is not approximately equal to {}",
            a.im,
            b.im
        );
    }};
    ($a:expr, $b:expr, $eps:expr) => {{
        let (a, b) = (&$a, &$b);
        let eps = $eps;
        assert!(
            ((a.re - b.re) as f64).abs() <= eps,
            "{} is not approximately equal to {}",
            a.re,
            b.re
        );
        assert!(
            ((a.im - b.im) as f64).abs() <= eps,
            "{} is not approximately equal to {}",
            a.im,
            b.im
        );
    }};
}
