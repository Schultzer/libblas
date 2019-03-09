use num_traits::{Float, NumAssignOps};

/// ROTMG
/// Note the version of ROTMG equal to [OpenBlas](https://github.com/xianyi/OpenBLAS/pull/1480)
/// based on an algorithm proposed by Tim Hopkins in 1998.
/// CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION [T] H WHICH ZEROS
/// THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(D1)*X1,SQRT(D2)*>    Y2)**T.
/// WITH PARAM(1)=FLAG, H HAS ONE OF THE FOLLOWING FORMS..
///
/// FLAG=-1.E0     FLAG=0.E0        FLAG=1.E0     FLAG=-2.E0
///
///   (H11  H12)    (1.E0  H12)    (H11  1.E0)    (1.E0  0.E0)
/// H=(          )    (          )    (          )    (          )
///   (H21  H22),   (H21  1.E0),   (-1.E0 H22),   (0.E0  1.E0).
/// LOCATIONS 2-4 OF PARAM CONTAIN H11,H21,H12, AND H22
/// RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
/// VALUE OF PARAM(1) ARE NOT STORED IN PARAM.)
///
/// THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
/// INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
/// OF D1 AND D2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
/// This is SROTMG and DROTMG comined in one function.
/// This is [SROTMG](http://www.netlib.org/lapack/explore-html/dd/d48/srotmg_8f.html) and [DROTMG](http://www.netlib.org/lapack/explore-html/df/deb/drotmg_8f.html) comined in one function.
#[inline]
pub fn rotmg<T: Float + NumAssignOps + num_traits::cast::FromPrimitive>(
    dd1: &mut T,
    dd2: &mut T,
    dx1: &mut T,
    dy1: &mut T,
    dparam: &mut [T],
) {
    let gam = T::from_f32(4096.0).unwrap();
    let gamsq = T::from_f32(16_777_216.0).unwrap();
    let rgamsq = T::from_f32(5.960_464_5E-8).unwrap();
    let two = T::one() + T::one();
    let du;
    let dp1;
    let dp2;
    let dq2;
    let dq1;
    let mut dh11 = T::zero();
    let mut dh21 = T::zero();
    let mut dh12 = T::zero();
    let mut dh22 = T::zero();
    let mut dflag = -T::one();
    let dtemp;

    if dd2.is_zero() || dy1.is_zero() {
        dflag = -two;
        dparam[0] = dflag;
        return;
    }

    if *dd1 < T::zero() {
        // dflag = -T::one();
        // dh11 = T::zero();
        // dh12 = T::zero();
        // dh21 = T::zero();
        // dh22 = T::zero();

        *dd1 = T::zero();
        *dd2 = T::zero();
        *dx1 = T::zero();
    } else if (dd1.is_zero() || dx1.is_zero()) && *dd2 > T::zero() {
        dflag = T::one();
        dh12 = T::one();
        dh21 = -T::one();
        *dx1 = *dy1;
        dtemp = *dd1;
        *dd1 = *dd2;
        *dd2 = dtemp;
    } else {
        dp2 = *dd2 * *dy1;
        if dp2.is_zero() {
            dflag = -two;
            dparam[0] = dflag;
            return;
        }
        dp1 = *dd1 * *dx1;
        dq2 = dp2 * *dy1;
        dq1 = dp1 * *dx1;
        if dq1.abs() > dq2.abs() {
            // dflag = T::zero();
            dh11 = T::one();
            dh22 = T::one();
            dh21 = -*dy1 / *dx1;
            dh12 = dp2 / dp1;

            du = T::one() - dh12 * dh21;
            if du > T::zero() {
                dflag = T::zero();
                *dd1 /= du;
                *dd2 /= du;
                *dx1 *= du;
            } else {
                dflag = -T::one();

                dh11 = T::zero();
                dh12 = T::zero();
                dh21 = T::zero();
                dh22 = T::zero();

                *dd1 = T::zero();
                *dd2 = T::zero();
                *dx1 = T::zero();
            }
        } else if dq2 < T::zero() {
            dflag = -T::one();

            dh11 = T::zero();
            dh12 = T::zero();
            dh21 = T::zero();
            dh22 = T::zero();

            *dd1 = T::zero();
            *dd2 = T::zero();
            *dx1 = T::zero();
        } else {
            dflag = T::one();
            dh21 = -T::one();
            dh12 = T::one();

            dh11 = dp1 / dp2;
            dh22 = *dx1 / *dy1;
            du = T::one() + dh11 * dh22;
            dtemp = *dd2 / du;

            *dd2 = *dd1 / du;
            *dd1 = dtemp;
            *dx1 = *dy1 * du;
        }

        while *dd1 <= rgamsq && !dd1.is_zero() {
            dflag = -T::one();
            *dd1 *= gam * gam;
            *dx1 /= gam;
            dh11 /= gam;
            dh12 /= gam;
        }
        while dd1.abs() > gamsq {
            dflag = -T::one();
            *dd1 /= gam * gam;
            *dx1 *= gam;
            dh11 *= gam;
            dh12 *= gam;
        }

        while dd2.abs() <= rgamsq && !dd2.is_zero() {
            dflag = -T::one();
            *dd2 *= gam * gam;
            dh21 /= gam;
            dh22 /= gam;
        }
        while dd2.abs() > gamsq {
            dflag = -T::one();
            *dd2 /= gam * gam;
            dh21 *= gam;
            dh22 *= gam;
        }
    }

    if dflag < T::zero() {
        dparam[1] = dh11;
        dparam[2] = dh21;
        dparam[3] = dh12;
        dparam[4] = dh22;
    } else if dflag.is_zero() {
        dparam[2] = dh21;
        dparam[3] = dh12;
    } else {
        dparam[1] = dh11;
        dparam[4] = dh22;
    }

    dparam[0] = dflag;
    return;
}


#[test]
fn test_rotmg() {
    let mut d1 = -4.0;
    let mut d2 = 2.0;
    let mut x1 = 3.0;
    let mut y1 = 9.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
            0.00024414062500,
            9.1552734375000000E-005
        ]
    );

    // GONUM TEST
    let mut d1 = 1600000000.;
    let mut d2 = 800000000.;
    let mut x1 = 8.;
    let mut y1 = 7.;
    let mut param = vec![1., 4096., -4096., 1., 4096.];
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 68.96627824858757);
    assert_eq!(d2, 34.483139124293785);
    assert_eq!(x1, 45312.);
    assert_eq!(y1, 7.0);
    assert_eq!(param, [-1., 4096., -3584., 1792., 4096.]);

    let mut d1 = 1f64 / (2 << 23) as f64;
    let mut d2 = 1f64 / (2 << 24) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 0.81818181818181812);
    assert_eq!(d2, 0.40909090909090906);
    assert_eq!(x1, 8.9518229166666674E-004);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [
            -1.0,
            2.4414062500000000E-004,
            -0.00016276041666666666,
            0.00008138020833333333,
            2.4414062500000000E-004
        ]
    );

    let mut d1 = 2f64 / (2 << 23) as f64;
    let mut d2 = 1f64 / (2 << 23) as f64;
    let mut x1 = 3.0;
    let mut y1 = 2.0;
    let mut param = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
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
    rotmg(&mut d1, &mut d2, &mut x1, &mut y1, &mut param);
    assert_eq!(d1, 1.6363636363636362);
    assert_eq!(d2, 13726813.090909090);
    assert_eq!(x1, 15018.666666666668);
    assert_eq!(y1, 2.0);
    assert_eq!(
        param,
        [-1.0, 4096.0, -0.66666666666666663, 1365.3333333333333, 1.0]
    );
}


