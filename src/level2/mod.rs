use num_traits::{Float, NumAssignOps};
use std::cmp::{max, min};

pub mod complex;

fn multiply<T: Float + NumAssignOps>(
    left: *mut T,
    right: T,
    len: isize,
    mut index: isize,
    inc: isize,
) {
    if right.is_zero() {
        zero(left, len, index, inc);
    } else {
        let mut i = 0;
        while i < len {
            unsafe {
                *left.offset(index) *= right;
            }
            index += inc;
            i += 1;
        }
    }
}

fn zero<T: Float + NumAssignOps>(left: *mut T, len: isize, mut index: isize, inc: isize) {
    let mut i = 0;
    while i < len {
        unsafe {
            *left.offset(index) = T::zero();
        }
        index += inc;
        i += 1;
    }
}

/// GBMV  performs one of the [T]-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
/// where alpha and beta are scalars, x and y are vectors and A is an m by n band [T], with kl sub-diagonals and ku super-diagonals.
/// This is [SGBMV](http://www.netlib.org/lapack/explore-html/d6/d46/sgbmv_8f.html) and [DGBMV](http://www.netlib.org/lapack/explore-html/d2/d3f/dgbmv_8f.html) comined in one function
#[inline]
pub fn gbmv<T: Float + NumAssignOps>(
    trans: char,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: T,
    a: *const T,
    lda: isize,
    x: *const T,
    incx: isize,
    beta: T,
    y: *mut T,
    incy: isize,
) {
    let mut info = 0;
    if trans != 'c' && trans != 'C' && trans != 'n' && trans != 'N' && trans != 't' && trans != 'T'
    {
        info = 1;
    } else if m < 0 {
        info = 2;
    } else if n < 0 {
        info = 3;
    } else if kl < 0 {
        info = 4;
    } else if ku < 0 {
        info = 5;
    } else if lda < (kl + ku + 1) {
        info = 8;
    } else if incx == 0 {
        info = 10;
    } else if incy == 0 {
        info = 13;
    }

    if info != 0 {
        panic!("gbmv {}", info);
    }

    if m == 0 || n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    let lenx = if trans == 'n' || trans == 'N' { n } else { m };
    let leny = if trans == 'n' || trans == 'N' { m } else { n };

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = ((-lenx) * incx) + incx
    };
    if incy < 0 {
        ky = ((-leny) * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, leny, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }

    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut iy = ky;
            // FIXME
            let k = ku - j; // apparently we want this to be negative sometimes.
            let aj = j * lda;
            let mut i = max(0, j - ku);
            while i < min(m, j + kl + 1) {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + k + i);
                }
                iy += incy;
                i += 1;
            }
            jx += incx;
            if j > ku - 1 {
                ky += incy
            };
            j += 1;
            //
        }
    } else {
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let mut tmp = T::zero();
            let mut ix = kx;
            let k = ku - j; // apparently we want this to be negative sometimes.
            let aj = j * lda;
            let mut i = max(0, j - ku);
            while i < min(m, j + kl + 1) {
                unsafe {
                    tmp += *a.offset(aj + k + i) * *x.offset(ix);
                }
                ix += incx;
                i += 1;
            }
            unsafe {
                *y.offset(jy) += alpha * tmp;
            }
            jy += incy;
            if j > ku - 1 {
                kx += incx
            };
            j += 1;
        }
    }
}

/// GEMV  performs one of the matrix-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
/// where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
/// This is [SGEMV](http://www.netlib.org/lapack/explore-html/db/d58/sgemv_8f.html) and [DGEMV](http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html) comined in one function
#[inline]
pub fn gemv<T: Float + NumAssignOps>(
    trans: char,
    m: isize,
    n: isize,
    alpha: T,
    a: *const T,
    lda: isize,
    x: *const T,
    incx: isize,
    beta: T,
    y: *mut T,
    incy: isize,
) {
    let mut info = 0;
    if trans != 'c' && trans != 'C' && trans != 'n' && trans != 'N' && trans != 't' && trans != 'T'
    {
        info = 1;
    } else if m < 0 {
        info = 2;
    } else if n < 0 {
        info = 3;
    } else if lda < max(1, m) {
        info = 6;
    } else if incx == 0 {
        info = 8;
    } else if incy == 0 {
        info = 11;
    }
    if info != 0 {
        panic!("gemv {}", info);
    }

    if m == 0 || n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    let lenx = if trans == 'n' || trans == 'N' { n } else { m };
    let leny = if trans == 'n' || trans == 'N' { m } else { n };

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-lenx * incx) + incx
    };
    if incy < 0 {
        ky = (-leny * incx) + incy
    };

    if !beta.is_one() {
        if incy == 1 && beta.is_zero() {
            zero(y, leny, ky, incy);
        } else {
            multiply(y, beta, leny, ky, incy)
        }
    }
    if alpha.is_zero() {
        return;
    }
    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + i);
                }
                iy += incy;
                i += 1;
            }
            jx += incx;
            j += 1;
        }
    } else {
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let mut tmp = T::zero();
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                unsafe {
                    tmp += *a.offset(aj + i) * *x.offset(ix);
                }
                ix += incx;
                i += 1;
            }
            unsafe {
                *y.offset(jy) += alpha * tmp;
            }
            jy += incy;
            j += 1;
        }
    }
}

/// SGER   performs the rank 1 operation
/// A := alpha*x*y**T + A,
/// where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix.
/// This is [SGER](http://www.netlib.org/lapack/explore-html/db/d5c/sger_8f.html) and [DGER](http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html) comined in one function
#[inline]
pub fn ger<T: Float + NumAssignOps>(
    m: isize,
    n: isize,
    alpha: T,
    x: *const T,
    incx: isize,
    y: *const T,
    incy: isize,
    a: *mut T,
    lda: isize,
) {
    let mut info = 0;
    if m < 0 {
        info = 1;
    } else if n < 0 {
        info = 2;
    } else if incx == 0 {
        info = 5;
    } else if incy == 0 {
        info = 7;
    } else if lda < max(1, m) {
        info = 9;
    }
    if info != 0 {
        panic!("ger {}", info);
    }

    if m == 0 || n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    let mut jy = 0;
    if incx < 0 {
        kx = (-m * incx) + incx
    };
    if incy < 0 {
        jy = (-n * incy) + incy
    };
    let mut j = 0;
    while j < n {
        let mut tmp;
        unsafe {
            tmp = *y.offset(jy);
        }
        if !tmp.is_zero() {
            tmp *= alpha;
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                unsafe {
                    *a.offset(aj + i) += *x.offset(ix) * tmp;
                }
                ix += incx;
                i += 1;
            }
        }
        jy += incy;
        j += 1;
    }
}

/// SBMV  performs the matrix-vector  operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n symmetric band matrix, with k super-diagonals.
/// This is [SSBMV](http://www.netlib.org/lapack/explore-html/d3/da1/ssbmv_8f.html) and [DSBMV](hhttp://www.netlib.org/lapack/explore-html/d8/d1e/dsbmv_8f.html) comined in one function
#[inline]
pub fn sbmv<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    k: isize,
    alpha: T,
    a: *const T,
    lda: isize,
    x: *const T,
    incx: isize,
    beta: T,
    y: *mut T,
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
    } else if k < 0 {
        info = 3;
    } else if lda < k + 1 {
        info = 6;
    } else if incx == 0 {
        info = 8;
    } else if incy == 0 {
        info = 11;
    }
    if info != 0 {
        panic!("sbmv {}", info);
    }

    if n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if incy < 0 {
        ky = (-n * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + k - j + i);
                    tmp2 += *a.offset(aj + k - j + i) * *x.offset(ix);
                }
                ix += incx;
                iy += incy;
                i += 1;
            }
            unsafe {
                *y.offset(jy) += tmp * *a.offset(aj + k) + alpha * tmp2;
            }
            jx += incx;
            jy += incy;
            j += 1;
            if j > k {
                kx += incx;
                ky += incy;
            }
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let tmp = unsafe { alpha * x.offset(jx).read() };
            let mut tmp2 = T::zero();
            unsafe {
                *y.offset(jy) += tmp * *a.offset(aj);
            }
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj - j + i);
                    tmp2 += *a.offset(aj - j + i) * *x.offset(ix);
                }
                i += 1;
            }
            unsafe {
                *y.offset(jy) += alpha * tmp2;
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// SPMV  performs the matrix-vector operation
/// y := alpha*A*x + beta*y,
// where alpha and beta are scalars, x and y are n element vectors and A is an n by n symmetric matrix, supplied in packed form.
/// This is [SSPMV](http://www.netlib.org/lapack/explore-html/d8/d68/sspmv_8f.html) and [DSPMV](http://www.netlib.org/lapack/explore-html/d4/d85/dspmv_8f.html) comined in one function
#[inline]
pub fn spmv<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    ap: *const T,
    x: *const T,
    incx: isize,
    beta: T,
    y: *mut T,
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2
    } else if incx == 0 {
        info = 6;
    } else if incy == 0 {
        info = 9;
    }
    if info != 0 {
        panic!("spmv {}", info);
    }

    if n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if incy < 0 {
        ky = (-n * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp;
            unsafe { tmp = alpha * *x.offset(jx) };
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j; //  kk += j - 1
            while k < kk - 1 {
                unsafe {
                    *y.offset(iy) += tmp * *ap.offset(k);
                    tmp2 += *ap.offset(k) * *x.offset(ix);
                }
                ix += incx;
                iy += incy;
                k += 1;
            }
            unsafe {
                *y.offset(jy) += tmp * *ap.offset(kk - 1) + alpha * tmp2;
            }
            jx += incx;
            jy += incy;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp;
            unsafe { tmp = alpha * *x.offset(jx) };
            let mut tmp2 = T::zero();
            unsafe { *y.offset(jy) += tmp * *ap.offset(kk) };
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk;
            kk += n - j;
            while k < kk {
                k += 1;
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *ap.offset(k);
                    tmp2 += *ap.offset(k) * *x.offset(ix);
                }
            }
            unsafe {
                *y.offset(jy) += alpha * tmp2;
            }
            jx += incx;
            jy += incy;
            kk += 1;
        }
    }
}

/// SPR    performs the symmetric rank 1 operation
/// A := alpha*x*x**T + A,
/// where alpha is a real scalar, x is an n element vector and A is an n by n symmetric matrix, supplied in packed form.
/// This is [SSPR](http://www.netlib.org/lapack/explore-html/d2/d9b/sspr_8f.html) and [DSPR](http://www.netlib.org/lapack/explore-html/dd/dba/dspr_8f.html) comined in one function
#[inline]
pub fn spr<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    x: *const T,
    incx: isize,
    ap: *mut T,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if incx == 0 {
        info = 6;
    }
    if info != 0 {
        panic!("spr {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    let mut kk = 0;
    let mut jx = kx;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut k = kk;
            kk += j;
            if !tmp.is_zero() {
                tmp *= alpha;
                let mut ix = kx;
                while k < kk {
                    unsafe { *ap.offset(k) += *x.offset(ix) * tmp };
                    ix += incx;
                    k += 1;
                }
            }
            jx += incx;
        }
    } else {
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut ix = jx;
            let mut k = kk;
            kk += n - j + 1;
            while k < kk {
                unsafe { *ap.offset(k) += *x.offset(ix) * tmp };
                ix += incx;
                k += 1;
            }
            jx += incx;
        }
    }
}

/// SPR2  performs the symmetric rank 2 operation
/// A := alpha*x*y**T + alpha*y*x**T + A,
/// where alpha is a scalar, x and y are n element vectors and A is an n by n symmetric matrix, supplied in packed form.
/// This is [SSPR2](http://www.netlib.org/lapack/explore-html/db/d3e/sspr2_8f.html) and [DSPR2](http://www.netlib.org/lapack/explore-html/dd/d9e/dspr2_8f.html) comined in one function
#[inline]
pub fn spr2<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    x: *const T,
    incx: isize,
    y: *const T,
    incy: isize,
    ap: *mut T,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if incx == 0 {
        info = 5;
    } else if incy == 0 {
        info = 7;
    }
    if info != 0 {
        panic!("spr2 {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if incy < 0 {
        ky = (-n * incy) + incy
    };

    let mut jx = kx;
    let mut jy = ky;
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp = unsafe { alpha * *y.offset(jy) };
            let tmp2 = unsafe { alpha * *x.offset(jx) };
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j;
            while k < kk {
                unsafe { *ap.offset(k) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                ix += incx;
                iy += incy;
                k += 1;
            }
            jx += incx;
            jy += incy;
        }
    } else {
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = unsafe { *y.offset(jy) };
            let mut tmp2 = unsafe { *x.offset(jx) };
            let mut k = kk;
            kk += n - j + 1;
            if !tmp.is_zero() || tmp2.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = jx;
                let mut iy = jy;
                while k < kk {
                    unsafe { *ap.offset(k) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                    ix += incx;
                    iy += incy;
                    k += 1;
                }
            }
            jx += incx;
            jy += incy;
        }
    }
}

/// SYMV  performs the matrix-vector  operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n symmetric matrix.
/// This is [SSYMV](http://www.netlib.org/lapack/explore-html/d2/d94/ssymv_8f.html) and [DSYMV](http://www.netlib.org/lapack/explore-html/d8/dbe/dsymv_8f.html) comined in one function
#[inline]
pub fn symv<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    a: *const T,
    lda: isize,
    x: *const T,
    incx: isize,
    beta: T,
    y: *mut T,
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if lda < max(1, n) {
        info = 5;
    } else if incx == 0 {
        info = 7;
    } else if incy == 0 {
        info = 9;
    }
    if info != 0 {
        panic!("symv {}", info);
    }

    if n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if incy < 0 {
        ky = (-n * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }
    let mut jx = kx;
    let mut jy = ky;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + i);
                    tmp2 += *a.offset(aj + i) * *x.offset(ix);
                }
                ix += incx;
                iy += incy;
                i += 1;
            }
            unsafe { *y.offset(jy) += tmp * *a.offset(aj + j) + alpha * tmp2 };
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2 = T::zero();
            let aj = j * lda;
            unsafe { *y.offset(jy) += tmp * *a.offset(aj + j) };
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + i);
                    tmp2 += *a.offset(aj + i) * *x.offset(ix)
                };
                i += 1;
            }
            unsafe { *y.offset(jy) += alpha * tmp2 };
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// SYR   performs the symmetric rank 1 operation
/// A := alpha*x*x**T + A,
/// where alpha is a real scalar, x is an n element vector and A is an n by n symmetric matrix.
/// This is [SSYR](http://www.netlib.org/lapack/explore-html/d6/dac/ssyr_8f.html) and [DSYR](http://www.netlib.org/lapack/explore-html/d3/d60/dsyr_8f.html) comined in one function
#[inline]
pub fn syr<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    x: *const T,
    incx: isize,
    a: *mut T,
    lda: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
    } else if incx == 0 {
        info = 5;
    } else if lda < max(1, n) {
        info = 7;
    }

    if info != 0 {
        panic!("syr {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    let mut jx = kx;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            let aj = j * lda;
            j += 1;
            let mut tmp = unsafe { *x.offset(jx) };
            if !tmp.is_zero() {
                tmp *= alpha;
                let mut ix = kx;
                let mut i = 0;
                while i < j {
                    // i < j + 1
                    unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp };
                    ix += incx;
                    i += 1;
                }
            }
            jx += incx;
        }
    } else {
        while j < n {
            let mut tmp = unsafe { *x.offset(jx) };
            if !tmp.is_zero() {
                tmp *= alpha;
                let aj = j * lda;
                let mut ix = jx;
                let mut i = j;
                while i < n {
                    unsafe {
                        let delta = *x.offset(ix) * tmp;
                        *a.offset(aj + i) += delta;
                    }
                    ix += incx;
                    i += 1;
                }
            }
            jx += incx;
            j += 1;
        }
    }
}

/// SYR2  performs the symmetric rank 2 operation
/// A := alpha*x*y**T + alpha*y*x**T + A,
/// where alpha is a scalar, x and y are n element vectors and A is an n by n symmetric matrix.
/// This is [SSYR2](http://www.netlib.org/lapack/explore-html/db/d99/ssyr2_8f.html) and [DSYR2](http://www.netlib.org/lapack/explore-html/db/d99/ssyr2_8f.html) comined in one function
#[inline]
pub fn syr2<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    x: *const T,
    incx: isize,
    y: *const T,
    incy: isize,
    a: *mut T,
    lda: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
    } else if incx == 0 {
        info = 3;
    } else if incy == 0 {
        info = 5;
    } else if lda < max(1, n) {
        info = 9;
    }

    if info != 0 {
        panic!("syr2 {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if incy < 0 {
        ky = (-n * incy) + incy
    };
    let mut jx = kx;
    let mut jy = ky;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            let mut tmp = unsafe { *y.offset(jy) };
            let mut tmp2 = unsafe { *x.offset(jx) };
            if !tmp2.is_zero() || !tmp.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = kx;
                let mut iy = ky;
                let aj = j * lda;
                let mut i = 0;
                while i < j + 1 {
                    unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                    ix += incx;
                    iy += incy;
                    i += 1;
                }
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let mut tmp = unsafe { *y.offset(jy) };
            let mut tmp2 = unsafe { *x.offset(jx) };
            if !tmp2.is_zero() || !tmp.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let aj = j * lda;
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j;
                while i < n {
                    unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                    ix += incx;
                    iy += incy;
                    i += 1;
                }
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// TBMV  performs one of the matrix-vector operations
/// x := A*x,   or   x := A**T*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular band matrix, with ( k + 1 ) diagonals.
/// This is [STMBV](http://www.netlib.org/lapack/explore-html/d6/d7d/stbmv_8f.html) and [DTMBV](http://www.netlib.org/lapack/explore-html/df/d29/dtbmv_8f.html) comined in one function
#[inline]
pub fn tbmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    k: isize,
    a: *const T,
    lda: isize,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if k < 0 {
        info = 5;
    } else if lda < k + 1 {
        info = 7;
    } else if incx == 0 {
        info = 9;
    }
    if info != 0 {
        panic!("tbmv {}", info);
    }

    if n == 0 {
        return;
    }
    let nounit = diag == 'n' || diag == 'N';
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = max(0, j - k);
                    while i < j {
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj + k - j + i) };
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        unsafe {
                            *x.offset(jx) *= *a.offset(aj + k);
                        }
                    }
                }
                jx += incx;
                j += 1;
                if j > k {
                    kx += incx;
                }
            }
        } else {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = min(n, j + k + 1);
                    while i > j + 1 {
                        i -= 1;
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj - j + i) };
                        ix -= incx;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) *= *a.offset(aj) };
                    }
                }
                jx -= incx;
                if (n - 1) - j >= k {
                    kx -= incx;
                }
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let aj = j * lda;
            let mut tmp = unsafe { *x.offset(jx) };
            kx -= incx;
            let mut ix = kx;
            if nounit {
                unsafe { tmp *= *a.offset(aj + k) };
            }
            let mut i = j;
            // FIXME figure out if there is a better way
            while i >= max(1, j - k + 1) {
                i -= 1;
                unsafe { tmp += *a.offset(aj + i + k - j) * *x.offset(ix) };
                ix -= incx;
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let mut tmp = unsafe { *x.offset(jx) };
            kx += incx;
            let mut ix = kx;
            if nounit {
                unsafe {
                    tmp *= *a.offset(aj);
                }
            }
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                unsafe { tmp += *a.offset(aj - j + i) * *x.offset(ix) };
                ix += incx;
                i += 1;
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            j += 1;
        }
    }
}

/// TBSV  solves one of the systems of equations
///  A*x = b,   or   A**T*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular band matrix, with ( k + 1 ) diagonals.
/// No test for singularity or near-singularity is included in this
/// routine. Such tests must be performed before calling this routine.
/// This is [STBSV](http://www.netlib.org/lapack/explore-html/d0/d1f/stbsv_8f.html) and [DTBSV](http://www.netlib.org/lapack/explore-html/d4/dcf/dtbsv_8f.html) comined in one function
#[inline]
pub fn tbsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    k: isize,
    a: *const T,
    lda: isize,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if k < 0 {
        info = 5;
    } else if lda < k + 1 {
        info = 7;
    } else if incx == 0 {
        info = 9;
    }
    if info != 0 {
        panic!("tbsv {}", info);
    }

    if n == 0 {
        return;
    }
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                kx -= incx;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let mut ix = kx;
                    let aj = j * lda;
                    if diag == 'n' || diag == 'N' {
                        unsafe {
                            *x.offset(jx) /= *a.offset(aj + k);
                        }
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut i = j;
                    while i > j - k {
                        i -= 1;
                        unsafe { *x.offset(ix) -= tmp * *a.offset(aj + k - j + i) };
                        ix -= incx;
                    }
                }
                jx -= incx;
            }
        } else {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                kx += incx;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let mut ix = kx;
                    let aj = j * lda;
                    if diag == 'n' || diag == 'N' {
                        unsafe { *x.offset(jx) /= *a.offset(aj) };
                    }
                    let mut i = j + 1;
                    let tmp = unsafe { *x.offset(jx) };
                    while i < min(n, j + k + 1) {
                        unsafe { *x.offset(ix) -= tmp * *a.offset(aj - j + i) };
                        ix += incx;
                        i += 1;
                    }
                }
                jx += incx;
                j += 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let aj = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                unsafe { tmp -= *a.offset(aj + k - j + i) * *x.offset(ix) };
                ix += incx;
                i += 1;
            }
            if diag == 'n' || diag == 'N' {
                unsafe {
                    tmp /= *a.offset(aj + k);
                }
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            j += 1;
            if j > k {
                kx += incx;
            }
        }
    } else {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let aj = j * lda;
            //FIXME Maybe we could do this in a diffrent way
            let mut i = min(n - 1, j + k);
            while i > j {
                unsafe { tmp -= *a.offset(aj - j + i) * *x.offset(ix) };
                ix -= incx;
                i -= 1;
            }
            //FIXME
            if diag == 'n' || diag == 'N' {
                unsafe {
                    tmp /= *a.offset(aj);
                }
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
            if (n - 1) - j >= k {
                kx -= incx;
            }
        }
    }
}

/// TPMV  performs one of the matrix-vector operations
/// x := A*x,   or   x := A**T*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// This is [STPMV](http://www.netlib.org/lapack/explore-html/db/db1/stpmv_8f.html) and [DTPMV](http://www.netlib.org/lapack/explore-html/dc/dcd/dtpmv_8f.html) comined in one function
#[inline]
pub fn tpmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    ap: *const T,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if incx == 0 {
        info = 7;
    }
    if info != 0 {
        panic!("tmpv {}", info);
    }

    if n == 0 {
        return;
    }
    let nounit = diag == 'n' || diag == 'N';
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut kk = 0;
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                j += 1;
                let mut k = kk;
                kk += j; // kk += j - 1;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = kx;
                    while k < kk - 1 {
                        unsafe {
                            *x.offset(ix) += tmp * *ap.offset(k);
                        }
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        unsafe {
                            *x.offset(jx) *= *ap.offset(kk - 1);
                        }
                    }
                }
                jx += incx;
            }
        } else {
            let mut kk = n * (n + 1) / 2;
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                //FIXME figure out a better way
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = kx;
                    let mut k = kk - 1;
                    while k > kk - (n - j) {
                        unsafe {
                            *x.offset(ix) += tmp * *ap.offset(k);
                        }
                        ix -= incx;
                        k -= 1;
                    }
                    if nounit {
                        unsafe {
                            *x.offset(jx) *= *ap.offset(kk - n + j);
                        }
                    }
                }
                jx -= incx;
                kk -= n - j;
                //
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut kk = n * (n + 1) / 2;
        let mut jx = kx + (n - 1) * incx;
        let mut j = n;
        while j >= 1 {
            kk -= j;
            j -= 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = jx;
            let mut k = kk + j;
            if nounit {
                unsafe {
                    tmp *= *ap.offset(k);
                }
            }
            while k > kk {
                ix -= incx;
                k -= 1;
                unsafe { tmp += *ap.offset(k) * *x.offset(ix) };
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
        }
    } else {
        let mut kk = 0;
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = jx;
            let mut k = kk;
            kk += n - j;
            if nounit {
                unsafe {
                    tmp *= *ap.offset(k);
                }
            }
            while k < kk {
                ix += incx;
                k += 1;
                unsafe { tmp += *ap.offset(k) * *x.offset(ix) };
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            kk += 1;
        }
    }
}

/// TPSV  solves one of the systems of equations
/// A*x = b,   or   A**T*x = b,
/// where b and x are n element vectors and A is an n by n unit, or routine. Such tests must be performed before calling this routine.
/// This is [STPSV](http://www.netlib.org/lapack/explore-html/d0/d7c/stpsv_8f.html) and [DTPSV](http://www.netlib.org/lapack/explore-html/d9/d84/dtpsv_8f.html) comined in one function
#[inline]
pub fn tpsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    ap: *const T,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4
    } else if incx == 0 {
        info = 7;
    }
    if info != 0 {
        panic!("tpsv {}", info);
    }

    if n == 0 {
        return;
    }
    let nounit = diag == 'n' || diag == 'N';
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut kk = n * (n + 1) / 2;
            let mut jx = kx + (n - 1) * incx;
            let mut j = n;
            while j >= 1 {
                kk -= j;
                j -= 1;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let mut k = kk + j;
                    if nounit {
                        unsafe {
                            *x.offset(jx) /= *ap.offset(k);
                        }
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = jx;
                    while k > kk {
                        ix -= incx;
                        k -= 1;
                        unsafe { *x.offset(ix) -= tmp * *ap.offset(k) };
                    }
                }
                jx -= incx;
            }
        } else {
            let mut kk = 0;
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                j += 1;
                let mut k = kk;
                kk += n - j;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    if nounit {
                        unsafe {
                            *x.offset(jx) /= *ap.offset(k);
                        }
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = jx;
                    while k < kk {
                        k += 1;
                        ix += incx;
                        unsafe { *x.offset(ix) -= tmp * *ap.offset(k) };
                    }
                }
                jx += incx;
                kk += 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut kk = 0;
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let mut k = kk;
            while k < kk + j - 1 {
                unsafe { tmp -= *ap.offset(k) * *x.offset(ix) };
                ix += incx;
                k += 1;
            }
            if nounit {
                // NOTE k == kk + j - 1
                unsafe {
                    tmp /= *ap.offset(k);
                }
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            kk += j;
        }
    } else {
        let mut kk = n * (n + 1) / 2;
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let mut k = kk - 1;
            kk -= n - j;
            while k > kk {
                unsafe { tmp -= *ap.offset(k) * *x.offset(ix) };
                ix -= incx;
                k -= 1;
            }
            if nounit {
                // NOTE k == kk - n + j
                unsafe {
                    tmp /= *ap.offset(k);
                }
            }
            unsafe {
                *x.offset(jx) = tmp;
            }
            jx -= incx;
        }
    }
}

/// TRMV  performs one of the matrix-vector operations
/// x := A*x,   or   x := A**T*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix.
/// This is [STRMV](http://www.netlib.org/lapack/explore-html/dc/d7e/dtrmv_8f.html) and [DTRMV](http://www.netlib.org/lapack/explore-html/de/d45/strmv_8f.html) comined in one function
#[inline]
pub fn trmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    a: *const T,
    lda: isize,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if lda < max(1, n) {
        info = 6;
    } else if incx == 0 {
        info = 8;
    }
    if info != 0 {
        panic!("trmv {}", info);
    }

    if n == 0 {
        return;
    }
    let nounit = diag == 'n' || diag == 'N';
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = 0;
                    while i < j {
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj + i) };
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) *= *a.offset(aj + j) };
                    }
                }
                jx += incx;
                j += 1;
            }
        } else {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let aj = j * lda;
                    let mut ix = kx;
                    let mut i = n;
                    while i > j + 1 {
                        i -= 1;
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj + i) };
                        ix -= incx;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) *= *a.offset(aj + j) };
                    }
                }
                jx -= incx;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut jx = kx + (n - 1) * incx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = jx;
            let aj = j * lda;
            let mut i = j;
            if nounit {
                unsafe {
                    tmp *= *a.offset(aj + j);
                }
            }
            while i >= 1 {
                i -= 1;
                ix -= incx;
                unsafe { tmp += *a.offset(aj + i) * *x.offset(ix) };
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = jx;
            let aj = j * lda;
            if nounit {
                unsafe {
                    tmp *= *a.offset(aj + j);
                }
            }
            let mut i = j + 1;
            while i < n {
                ix += incx;
                unsafe { tmp += *a.offset(aj + i) * *x.offset(ix) };
                i += 1;
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            j += 1;
        }
    }
}

/// TRSV  solves one of the systems of equations
/// A*x = b,   or   A**T*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular matrix.
/// No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.
/// This is [STRSV](http://www.netlib.org/lapack/explore-html/d0/d2a/strsv_8f.html) and [DTRSV](http://www.netlib.org/lapack/explore-html/d6/d96/dtrsv_8f.html) comined in one function
#[inline]
pub fn trsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    a: *const T,
    lda: isize,
    x: *mut T,
    incx: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 2;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if lda < max(1, n) {
        info = 6;
    } else if incx == 0 {
        info = 8;
    }
    if info != 0 {
        panic!("trsv {}", info);
    }

    if n == 0 {
        return;
    }
    let nounit = diag == 'n' || diag == 'N';
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx + (n - 1) * incx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let aj = j * lda;
                    if nounit {
                        unsafe { *x.offset(jx) /= *a.offset(aj + j) };
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = jx;
                    let mut i = j;
                    while i >= 1 {
                        i -= 1;
                        ix -= incx;
                        unsafe { *x.offset(ix) -= tmp * *a.offset(aj + i) };
                    }
                }
                jx -= incx;
            }
        } else {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let aj = j * lda;
                    if nounit {
                        unsafe { *x.offset(jx) /= *a.offset(aj + j) };
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = jx;
                    let mut i = j + 1;
                    while i < n {
                        ix += incx;
                        unsafe {
                            *x.offset(ix) -= tmp * *a.offset(aj + i);
                        }
                        i += 1;
                    }
                }
                jx += incx;
                j += 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                unsafe {
                    tmp -= *a.offset(aj + i) * *x.offset(ix);
                }
                ix += incx;
                i += 1;
            }
            if nounit {
                unsafe {
                    tmp /= *a.offset(aj + i);
                }
            }
            unsafe {
                *x.offset(jx) = tmp;
            }
            jx += incx;
            j += 1;
        }
    } else {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            let aj = j * lda;
            let mut i = n;
            while i > j + 1 {
                i -= 1;
                unsafe {
                    tmp -= *a.offset(aj + i) * *x.offset(ix);
                }
                ix -= incx;
            }
            if nounit {
                unsafe {
                    tmp /= *a.offset(aj + j);
                }
            }
            unsafe {
                *x.offset(jx) = tmp;
            }
            jx -= incx;
        }
    }
}
