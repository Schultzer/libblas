// use std::arch::x86_64::*;
use num_traits::{Float, NumAssignOps};

pub mod complex;

/// IAMAX finds the index of the first element having maximum absolute value.
/// This is [ISAMAX](http://www.netlib.org/lapack/explore-html/d6/d44/isamax_8f.html) and [IDAMAX](http://www.netlib.org/lapack/explore-html/dd/de0/idamax_8f.html) comined in one function
#[inline]
pub fn iamax<T: Float + NumAssignOps>(n: usize, x: &[T], incx: usize) -> usize {
    let mut iamax = 0;
    if n < 1 || incx <= 0 {
        return iamax;
    }

    iamax = 1;
    if n == 1 {
        return iamax;
    }
    let mut i = 2;
    let mut max = x[0].abs();
    if incx == 1 {
        while i < n {
            let tmp = x[i].abs();
            i += 1;
            if tmp > max {
                iamax = i;
                max = tmp;
            }
        }
    } else {
        let mut ix = 1 + incx;
        while i < n {
            let tmp = x[ix].abs();
            ix += incx;
            i += 1;
            if tmp > max {
                iamax = i;
                max = tmp;
            }
        }
    }
    iamax
}

/// ASUM takes the sum of the absolute values.
/// This is [SASUM](http://www.netlib.org/lapack/explore-html/df/d1f/sasum_8f.html) and [DASUM](http://www.netlib.org/lapack/explore-html/de/d05/dasum_8f.html) comined in one function
#[inline]
pub fn asum<T: Float + NumAssignOps>(n: usize, x: &[T], incx: usize) -> T {
    let mut asum = T::zero();
    if n == 0 || incx == 0 {
        return asum;
    }
    let mut i = 0;
    let nincx = n * incx;
    while i < nincx {
        asum += x[i].abs();
        i += incx;
    }
    asum
}

/// AXPY constant times a vector plus a vector.
/// This is [SAXPY](http://www.netlib.org/lapack/explore-html/d8/daf/saxpy_8f.html) and [DAXPY](http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html) comined in one function
#[inline]
pub fn axpy<T: Float + NumAssignOps>(
    n: usize,
    a: T,
    x: &[T],
    incx: isize,
    y: &mut [T],
    incy: isize,
) {
    if n == 0 || a.is_zero() {
        return;
    };
    let mut i = 0;
    if incx == 1 && incy == 1 {
        let m = n % 4;
        if m != 0 {
            while i < m {
                y[i] += a * x[i];
                i += 1;
            }
        }
        if n < 4 {
            return;
        };
        let mut i = m;
        while i < n {
            y[i] += a * x[i];
            i += 1;
            y[i] += a * x[i];
            i += 1;
            y[i] += a * x[i];
            i += 1;
            y[i] += a * x[i];
            i += 1;
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 {
            ix = (-(n as isize) * incx) + incx;
        }
        if incy < 0 {
            iy = (-(n as isize) * incy) + incy;
        }
        while i < n {
            y[iy as usize] += a * x[ix as usize];
            ix += incx;
            iy += incy;
            i += 1;
        }
    }
}

/// NRM2 returns the euclidean norm of a vector via the function name, so that NRM2 := sqrt( x'*x ).
/// This is [SNRM2](http://www.netlib.org/lapack/explore-html/d7/df1/snrm2_8f.html) and [DNRM2](http://www.netlib.org/lapack/explore-html/da/d7f/dnrm2_8f.html) comined in one function
#[inline]
pub fn nrm2<T: Float + NumAssignOps>(n: usize, x: &[T], incx: usize) -> T {
    let mut scale = T::zero();
    let mut ssq = T::one();
    if n == 0 || incx == 0 {
        return scale;
    }
    if n == 1 {
        return x[0].abs();
    }
    //
    //  The following loop is equivalent to this call to the LAPACK auxiliary routine:
    //  CALL SLASSQ( N, X, INCX, SCALE, SSQ )
    //
    let mut ix = 0;
    while ix < (n - 1) * incx {
        let r = x[ix];
        if r != T::zero() {
            let tmp = r.abs();
            let ratio = scale / tmp;
            let ratiop2 = ratio * ratio;
            if scale < tmp {
                ssq = T::one() + ssq * ratiop2;
                scale = tmp;
            } else {
                ssq += T::one() / ratiop2;
            }
        }
        ix += incx
    }
    scale * ssq.sqrt()
}

/// COPY copies a vector, x, to a vector, y. uses unrolled loops for increments equal to 1.
/// This is [SCOPY](http://www.netlib.org/lapack/explore-html/de/dc0/scopy_8f.html) and [DCOPY](http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html) comined in one function
#[inline]
pub fn copy<T: Float + NumAssignOps>(n: usize, x: &[T], incx: isize, y: &mut [T], incy: isize) {
    if n == 0 {
        return;
    }
    let mut i = 0;
    if incx == 1 && incy == 1 {
        let m = n % 7;
        if m != 0 {
            while i < m {
                y[i] = x[i];
                i += 1;
            }
        }
        if n < 7 {
            return;
        }
        let mut i = m;
        while i < n {
            y[i] = x[i]; // 1
            i += 1;
            y[i] = x[i];
            i += 1;
            y[i] = x[i];
            i += 1;
            y[i] = x[i];
            i += 1;
            y[i] = x[i];
            i += 1;
            y[i] = x[i];
            i += 1;
            y[i] = x[i]; // 7
            i += 1;
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 {
            ix = (-(n as isize) * incx) + incx;
        }
        if incy < 0 {
            iy = (-(n as isize) * incy) + incy;
        }
        while i < n {
            y[iy as usize] = x[ix as usize];
            ix += incx;
            iy += incy;
            i += 1;
        }
    }
}

/// DOT forms the dot product of two vectors. uses unrolled loops for increments equal to one.
/// This is [SDOT](http://www.netlib.org/lapack/explore-html/d0/d16/sdot_8f.html) and [DDOT](http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html) comined in one function
#[inline]
pub fn dot<T: Float + NumAssignOps>(n: usize, x: &[T], incx: isize, y: &[T], incy: isize) -> T {
    let mut dot = T::zero();
    if n == 0 {
        return dot;
    }
    let mut i = 0;
    if incx == 1 && incy == 1 {
        let m = n % 5;
        if m != 0 {
            while i < m {
                dot += x[i] * y[i];
                i += 1;
            }
            if n < 5 {
                return dot;
            }
        }
        let mut i = m;
        while i < n {
            dot += x[i] * y[i];
            i += 1;
            dot += x[i] * y[i];
            i += 1;
            dot += x[i] * y[i];
            i += 1;
            dot += x[i] * y[i];
            i += 1;
            dot += x[i] * y[i];
            i += 1;
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 {
            ix = (-(n as isize) * incx) + incx;
        }
        if incy < 0 {
            iy = (-(n as isize) * incy) + incy;
        }
        while i < n {
            dot += x[ix as usize] * y[iy as usize];
            ix += incx;
            iy += incy;
            i += 1;
        }
    }
    dot
}

/// DDOT
/// Compute the inner product of two vectors with extended precision accumulation.
/// This is [SDSDOT](http://www.netlib.org/lapack/explore-html/d9/d47/sdsdot_8f.html) and [DSDOT](http://www.netlib.org/lapack/explore-html/dc/d01/dsdot_8f.html) comined in one function
#[inline]
pub fn ddot<T: Float + NumAssignOps>(
    n: usize,
    b: T,
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
) -> f64 {
    let mut dot = b.to_f64().unwrap();
    if n == 0 {
        return dot;
    }
    //
    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-(n as isize) * incx) + incx;
    }
    if incy < 0 {
        iy = (-(n as isize) * incy) + incy;
    }
    while i < n {
        dot += (x[ix as usize] * y[iy as usize]).to_f64().unwrap();
        ix += incx;
        iy += incy;
        i += 1;
    }
    dot
}

/// ROT applies a plane rotation
/// This is [SROT](http://www.netlib.org/lapack/explore-html/db/d6c/srot_8f.html) and [DROT](http://www.netlib.org/lapack/explore-html/dc/d23/drot_8f.html) comined in one function
#[inline]
pub fn rot<T: Float + NumAssignOps>(
    n: usize,
    x: &mut [T],
    incx: isize,
    y: &mut [T],
    incy: isize,
    c: T,
    s: T,
) {
    if n == 0 {
        return;
    }
    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-(n as isize) * incx) + incx;
    }
    if incy < 0 {
        iy = (-(n as isize) * incy) + incy;
    }
    while i < n {
        let rot = c * x[ix as usize] + s * y[iy as usize];
        y[iy as usize] = -s * x[ix as usize] + c * y[iy as usize];
        x[ix as usize] = rot;
        i += 1;
        ix += incx;
        iy += incy;
    }
}

/// ROTG construct givens plane rotation.
/// This is [SROTG](http://www.netlib.org/lapack/explore-html/d7/d26/srotg_8f.html) and [DROTG](http://www.netlib.org/lapack/explore-html/de/d13/drotg_8f.html) comined in one function.
#[inline]
pub fn rotg<T: Float + NumAssignOps>(a: &mut T, b: &mut T, c: &mut T, s: &mut T) {
    let roe = if a.abs() > b.abs() { *a } else { *b };
    let scale = a.abs() + b.abs();
    if scale == T::zero() {
        *c = T::one();
        *s = T::zero();
        *a = T::zero();
        *b = T::zero();
        return;
    }
    let v1 = *a / scale;
    let v2 = *b / scale;
    let mut r = scale * (v1 * v1 + v2 * v2).sqrt();
    r = if roe >= T::zero() { r } else { -r };
    *c = *a / r;
    *s = *b / r;
    let mut z = T::one();
    if a.abs() > b.abs() {
        z = *s
    }
    if a.abs() <= b.abs() && *c != T::zero() {
        z = T::one() / *c
    }
    *a = r;
    *b = z;
}

/// ROTM
/// APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N [T]
/// This is [SROTM](http://www.netlib.org/lapack/explore-html/d6/d0f/srotm_8f.html) and [DROTM](http://www.netlib.org/lapack/explore-html/d8/d7b/drotm_8f.html) comined in one function.
#[inline]
pub fn rotm<T: Float + NumAssignOps>(
    n: usize,
    x: &mut [T],
    incx: isize,
    y: &mut [T],
    incy: isize,
    param: &[T],
) {
    let flag = param[0];
    let mut i = 0;
    let h11 = param[1];
    let h12 = param[3];
    let h21 = param[2];
    let h22 = param[4];
    let two = T::one() + T::one();
    if n == 0 || flag + two == T::zero() {
        return;
    };
    if incx == incy && incx > 0 {
        let nsteps = n * incx as usize;
        if flag < T::zero() {
            while i < nsteps {
                let w = x[i];
                let z = y[i];

                x[i] = w * h11 + z * h12;
                y[i] = w * h21 + z * h22;
                i += incx as usize;
            }
        } else if flag.is_zero() {
            while i < nsteps {
                let w = x[i];
                let z = y[i];
                x[i] = w + z * h12;
                y[i] = w * h21 + z;
                i += incx as usize;
            }
        } else {
            while i < nsteps {
                let w = x[i];
                let z = y[i];
                x[i] = w * h11 + z;
                y[i] = -w + h22 * z;
                i += incx as usize;
            }
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 {
            ix = (-(n as isize) * incx) + incx;
        }
        if incy < 0 {
            iy = (-(n as isize) * incy) + incy;
        }
        if flag < T::zero() {
            while i < n {
                let w = x[ix as usize];
                let z = y[iy as usize];

                x[ix as usize] = w * h11 + z * h12;
                y[iy as usize] = w * h21 + z * h22;
                ix += incx;
                iy += incy;
                i += 1;
            }
        } else if flag.is_zero() {
            while i < n {
                let w = x[ix as usize];
                let z = y[iy as usize];
                x[ix as usize] = w + z * h12;
                y[iy as usize] = w * h21 + z;
                ix += incx;
                iy += incy;
                i += 1;
            }
        } else {
            while i < n {
                let w = x[ix as usize];
                let z = y[iy as usize];
                x[ix as usize] = w * h11 + z;
                y[iy as usize] = -w + h22 * z;
                ix += incx;
                iy += incy;
                i += 1;
            }
        }
    }
}

/// ROTMG
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

/// NETLIB VERSERION;
// pub fn rotmg<T: Float + NumAssignOps + num_traits::cast::FromPrimitive>(
//     d1: &mut T,
//     d2: &mut T,
//     x1: &mut T,
//     y1: &mut T,
//     param: &mut [T],
// ) {
//     // https://github.com/xianyi/OpenBLAS/pull/1480
//     let mut flag = param[0];
//     let mut h11 = param[1];
//     let mut h12 = param[3];
//     let mut h21 = param[2];
//     let mut h22 = param[4];
//     let gam = T::from_f32(4096.0).unwrap();
//     let gamsq = T::from_f32(16_777_216.0).unwrap();
//     let rgamsq = T::from_f32(5.960_464_5E-8).unwrap();

//     if *d1 < T::zero() {
//         flag = -T::one();
//         h11 = T::zero();
//         h12 = T::zero();
//         h21 = T::zero();
//         h22 = T::zero();
//         *d1 = T::zero();
//         *d2 = T::zero();
//         *x1 = T::zero();
//     } else {
//         let p2 = *d2 * *y1;
//         if p2 == T::zero() {
//             param[0] = -(T::one() + T::one());
//             return;
//         }
//         let p1 = *d1 * *x1;
//         let q2 = p2 * *y1;
//         let q1 = p1 * *x1;
//         if q1.abs() > q2.abs() {
//             h21 = -*y1 / *x1;
//             h12 = p2 / p1;
//             let u = T::one() - h12 * h21;
//             if u > T::zero() {
//                 flag = T::zero();
//                 *d1 /= u;
//                 *d2 /= u;
//                 *x1 *= u;
//             }
//         } else if q2 < T::zero() {
//             flag = -T::one();
//             h11 = T::zero();
//             h12 = T::zero();
//             h21 = T::zero();
//             h22 = T::zero();
//             *d1 = T::zero();
//             *d2 = T::zero();
//             *x1 = T::zero();
//         } else {
//             flag = T::one();
//             h11 = p1 / p2;
//             h22 = *x1 / *y1;
//             let u = T::one() + h11 * h22;
//             let tmp = *d2 / u;
//             *d2 = *d1 / u;
//             *d1 = tmp;
//             *x1 = *y1 * u;
//         }
//         if *d1 != T::zero() {
//             while (*d1 <= rgamsq) || (*d1 >= gamsq) {
//                 if flag.is_zero() {
//                     h11 = T::one();
//                     h22 = T::one();
//                     flag = -T::one();
//                 } else {
//                     h21 = -T::one();
//                     h12 = T::one();
//                     flag = -T::one();
//                 }
//                 if *d1 <= rgamsq {
//                     *d1 *= gamsq;
//                     *x1 /= gam;
//                     h11 /= gam;
//                     h12 /= gam;
//                 } else {
//                     *d1 /= gamsq;
//                     *x1 *= gam;
//                     h11 *= gam;
//                     h12 *= gam;
//                 }
//             }
//         }
//         if *d2 != T::zero() {
//             while (d2.abs() <= rgamsq) || (d2.abs() >= gamsq) {
//                 if flag.is_zero() {
//                     h11 = T::one();
//                     h22 = T::one();
//                     flag = -T::one();
//                 } else {
//                     h21 = -T::one();
//                     h12 = T::one();
//                     flag = -T::one();
//                 }
//                 if d2.abs() <= rgamsq {
//                     *d2 *= gamsq;
//                     h21 /= gam;
//                     h22 /= gam;
//                 } else {
//                     *d2 /= gamsq;
//                     h21 *= gam;
//                     h22 *= gam;
//                 }
//             }
//         }
//     }
//     if flag < T::zero() {
//         param[1] = h11;
//         param[2] = h21;
//         param[3] = h12;
//         param[4] = h22;
//     } else if flag.is_zero() {
//         param[2] = h21;
//         param[3] = h12;
//     } else {
//         param[1] = h11;
//         param[4] = h22;
//     }
//     param[0] = flag;
// }

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

/// SCAL scales a vector by a constant. uses unrolled loops for increment equal to 1.
/// This is [SSCAL](http://www.netlib.org/lapack/explore-html/d9/d04/sscal_8f.html) and [DSCAL](http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html) comined in one function.
#[inline]
pub fn scal<T: Float + NumAssignOps>(n: usize, a: T, x: &mut [T], incx: usize) {
    if n <= 0 || incx <= 0 {
        return;
    };
    if a.is_one() {
        return;
    };
    let mut i = 0;
    if incx == 1 {
        let m = n % 5;
        if m != 0 {
            while i < m {
                x[i] *= a;
                i += 1;
            }
            if n < 5 {
                return;
            };
        }
        let mut mp1 = m;
        while mp1 < n {
            x[mp1] *= a;
            mp1 += 1;
            x[mp1] *= a;
            mp1 += 1;
            x[mp1] *= a;
            mp1 += 1;
            x[mp1] *= a;
            mp1 += 1;
            x[mp1] *= a;
            mp1 += 1;
        }
    } else {
        let nincx = n * incx;
        while i < nincx {
            x[i] *= a;
            i += incx;
        }
    }
}

/// SWAP interchanges two vectors. uses unrolled loops for increment equal to 1.
/// This is [SSWAP](http://www.netlib.org/lapack/explore-html/d9/da9/sswap_8f.html) and [DSWAP](http://www.netlib.org/lapack/explore-html/db/dd4/dswap_8f.html) comined in one function.
#[inline]
pub fn swap<T: Float + NumAssignOps>(n: usize, x: &mut [T], incx: isize, y: &mut [T], incy: isize) {
    if n <= 0 {
        return;
    };
    if incx == 1 && incy == 1 {
        let m = n % 3;
        if m != 0 {
            let mut i = 0;
            while i < m {
                std::mem::swap(&mut x[i], &mut y[i]);
                i += 1;
            }
            if n < 3 {
                return;
            };
        }
        let mut mp1 = m;
        while mp1 < n {
            std::mem::swap(&mut x[mp1], &mut y[mp1]);
            mp1 += 1;
            std::mem::swap(&mut x[mp1], &mut y[mp1]);
            mp1 += 1;
            std::mem::swap(&mut x[mp1], &mut y[mp1]);
            mp1 += 1;
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        let mut i = 0;

        if incx < 0 {
            ix = (-(n as isize) * incx) + incx;
        }
        if incy < 0 {
            iy = (-(n as isize) * incy) + incy;
        }
        while i < n {
            std::mem::swap(&mut x[ix as usize], &mut y[iy as usize]);
            ix += incx;
            iy += incy;
            i += 1;
        }
    }
}
