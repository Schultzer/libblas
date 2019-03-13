// use std::arch::x86_64::*;
use num_traits::{Float, NumAssignOps};
use std::cmp::{max, min};

pub mod complex;

fn multiply<T: Float + NumAssignOps>(
    left: &mut [T],
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
            left[index as usize] *= right;
            index += inc;
            i += 1;
        }
    }
}

fn zero<T: Float + NumAssignOps>(left: &mut [T], len: isize, mut index: isize, inc: isize) {
    let mut i = 0;
    while i < len {
        left[index as usize] = T::zero();
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
    m: usize,
    n: usize,
    kl: usize,
    ku: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
    incy: isize,
) {
    let mut info = 0;
    if trans != 'c' && trans != 'C' && trans != 'n' && trans != 'N' && trans != 't' && trans != 'T'
    {
        info = 1;
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
        kx = (-(lenx as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(leny as isize) * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, leny as isize, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }

    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut iy = ky;
            // FIXME
            let k = ku as isize - j as isize; // apparently we want this to be negative sometimes.
            let aj = (j * lda) as isize;
            let mut i = j.saturating_sub(ku); // MAX(0, J-KU) using saturating_sub it will always be 0
            while i < min(m, j + kl + 1) {
                y[iy as usize] += tmp * a[(aj + k + i as isize) as usize];
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
            let k = ku as isize - j as isize; // apparently we want this to be negative sometimes.
            let aj = (j * lda) as isize;
            let mut i = j.saturating_sub(ku); // MAX(0, J-KU) using saturating_sub it will always be 0
            while i < min(m, j + kl + 1) {
                tmp += a[(aj + k + i as isize) as usize] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            y[jy as usize] += alpha * tmp;
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
    m: usize,
    n: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
    incy: isize,
) {
    let mut info = 0;
    if trans != 'c' && trans != 'C' && trans != 'n' && trans != 'N' && trans != 't' && trans != 'T'
    {
        info = 1;
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
        kx = (-(lenx as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(leny as isize) * incx) + incy
    };

    if !beta.is_one() {
        if incy == 1 && beta.is_zero() {
            zero(y, leny as isize, ky, incy);
        } else {
            multiply(y, beta, leny as isize, ky, incy)
        }
    }
    if alpha.is_zero() {
        return;
    }
    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                y[iy as usize] += tmp * a[aj + i];
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
                tmp += a[aj + i] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            y[jy as usize] += alpha * tmp;
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
    m: usize,
    n: usize,
    alpha: T,
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    a: &mut [T],
    lda: usize,
) {
    let mut info = 0;
    if incx == 0 {
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
        kx = (-(m as isize) * incx) + incx
    };
    if incy < 0 {
        jy = (-(n as isize) * incy) + incy
    };
    let mut j = 0;
    while j < n {
        let mut tmp = y[jy as usize];
        if !tmp.is_zero() {
            tmp *= alpha;
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                a[aj + i] += x[ix as usize] * tmp;
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
    n: usize,
    k: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
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
        kx = (-(n as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(n as isize) * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n as isize, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = j.saturating_sub(k); // MAX(0, j - k)
            while i < j {
                y[iy as usize] += tmp * a[aj + k - j + i];
                tmp2 += a[aj + k - j + i] * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[(aj + k) as usize] + alpha * tmp2;
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
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            y[jy as usize] += tmp * a[aj];
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[aj - j + i];
                tmp2 += a[aj - j + i] * x[ix as usize];
                i += 1;
            }
            y[jy as usize] += alpha * tmp2;
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
    n: usize,
    alpha: T,
    ap: &[T],
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
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
        kx = (-(n as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(n as isize) * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n as isize, ky, incy)
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
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j;
            while k < kk - 1 {
                y[iy as usize] += tmp * ap[k];
                tmp2 += ap[k] * x[ix as usize];
                ix += incx;
                iy += incy;
                k += 1;
            }
            y[jy as usize] += tmp * ap[kk - 1] + alpha * tmp2;
            jx += incx;
            jy += incy;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            y[jy as usize] += tmp * ap[kk];
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk;
            kk += n - j;
            while k < kk {
                k += 1;
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * ap[k];
                tmp2 += ap[k] * x[ix as usize];
            }
            y[jy as usize] += alpha * tmp2;
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
    n: usize,
    alpha: T,
    x: &[T],
    incx: isize,
    ap: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };
    let mut kk = 0;
    let mut jx = kx;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = x[jx as usize];
            let mut k = kk;
            kk += j;
            if !tmp.is_zero() {
                tmp *= alpha;
                let mut ix = kx;
                while k < kk {
                    ap[k] += x[ix as usize] * tmp;
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
            let tmp = alpha * x[jx as usize];
            let mut ix = jx;
            let mut k = kk;
            kk += n - j + 1;
            while k < kk {
                ap[k] += x[ix as usize] * tmp;
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
    n: usize,
    alpha: T,
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    ap: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(n as isize) * incy) + incy
    };

    let mut jx = kx;
    let mut jy = ky;
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp = alpha * y[jy as usize];
            let tmp2 = alpha * x[jx as usize];
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j;
            while k < kk {
                ap[k] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
            let mut tmp = y[jy as usize];
            let mut tmp2 = x[jx as usize];
            let mut k = kk;
            kk += n - j + 1;
            if !tmp.is_zero() || tmp2.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = jx;
                let mut iy = jy;
                while k < kk {
                    ap[k] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
    n: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(n as isize) * incy) + incy
    };

    if !beta.is_one() {
        multiply(y, beta, n as isize, ky, incy)
    }
    if alpha.is_zero() {
        return;
    }
    let mut jx = kx;
    let mut jy = ky;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                y[iy as usize] += tmp * a[aj + i];
                tmp2 += a[aj + i] * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[aj + j] + alpha * tmp2;
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let aj = j * lda;
            y[jy as usize] += tmp * a[aj + j];
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[aj + i];
                tmp2 += a[aj + i] * x[ix as usize];
                i += 1;
            }
            y[jy as usize] += alpha * tmp2;
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
    n: usize,
    alpha: T,
    x: &[T],
    incx: isize,
    a: &mut [T],
    lda: usize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
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
        kx = (-(n as isize) * incx) + incx
    };
    let mut jx = kx;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            // FIXME check the ASM else move back into if statement
            let aj = j * lda;
            j += 1;
            //
            let mut tmp = x[jx as usize];
            if !tmp.is_zero() {
                tmp *= alpha;
                let mut ix = kx;
                let mut i = 0;
                while i < j {
                    // i < j + 1
                    a[aj + i] += x[ix as usize] * tmp;
                    ix += incx;
                    i += 1;
                }
            }
            jx += incx;
        }
    } else {
        while j < n {
            let mut tmp = x[jx as usize];
            if !tmp.is_zero() {
                tmp *= alpha;
                let aj = j * lda;
                let mut ix = jx;
                let mut i = j;
                while i < n {
                    let delta = x[ix as usize] * tmp;
                    a[aj + i] += delta;
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
    n: usize,
    alpha: T,
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    a: &mut [T],
    lda: usize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
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
        kx = (-(n as isize) * incx) + incx
    };
    if incy < 0 {
        ky = (-(n as isize) * incy) + incy
    };
    let mut jx = kx;
    let mut jy = ky;
    let mut j = 0;
    if uplo == 'u' || uplo == 'U' {
        while j < n {
            let mut tmp = y[jy as usize];
            let mut tmp2 = x[jx as usize];
            if !tmp2.is_zero() || !tmp.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = kx;
                let mut iy = ky;
                let aj = j * lda;
                let mut i = 0;
                // FIXME leave this as is for now check the ASM in syr if there is anything to move j * lda out of the if statement
                while i < j + 1 {
                    a[aj + i] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
            let mut tmp = y[jy as usize];
            let mut tmp2 = x[jx as usize];
            if !tmp2.is_zero() || !tmp.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let aj = j * lda;
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j;
                while i < n {
                    a[aj + i] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
    n: usize,
    k: usize,
    a: &[T],
    lda: usize,
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = j.saturating_sub(k);
                    while i < j {
                        x[ix as usize] += tmp * a[aj + k - j + i];
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[aj + k];
                    }
                }
                jx += incx;
                j += 1;
                if j > k {
                    kx += incx;
                }
            }
        } else {
            kx += (n as isize - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = min(n, j + k + 1);
                    while i > j + 1 {
                        i -= 1;
                        x[ix as usize] += tmp * a[aj - j + i];
                        ix -= incx;
                    }
                    if nounit {
                        x[jx as usize] *= a[aj];
                    }
                }
                jx -= incx;
                if (n - 1) - j >= k {
                    kx -= incx;
                }
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        kx += (n as isize - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let aj = j * lda;
            let mut tmp = x[jx as usize];
            kx -= incx;
            let mut ix = kx;
            if nounit {
                tmp *= a[aj + k];
            }
            let mut i = j;
            // FIXME figure out if there is a better way
            while i >= max(1, j.saturating_sub(k) + 1) {
                i -= 1;
                tmp += a[aj + i + k - j] * x[ix as usize];
                ix -= incx;
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let mut tmp = x[jx as usize];
            kx += incx;
            let mut ix = kx;
            if nounit {
                tmp *= a[(aj) as usize];
            }
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                tmp += a[aj - j + i] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            x[jx as usize] = tmp;
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
    n: usize,
    k: usize,
    a: &[T],
    lda: usize,
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            kx += (n as isize - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                kx -= incx;
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    if diag == 'n' || diag == 'N' {
                        x[jx as usize] /= a[aj + k];
                    }
                    let tmp = x[jx as usize];
                    let mut i = j;
                    while i > j.saturating_sub(k) {
                        i -= 1;
                        x[ix as usize] -= tmp * a[aj + k - j + i];
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
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    if diag == 'n' || diag == 'N' {
                        x[jx as usize] /= a[aj];
                    }
                    let mut i = j + 1;
                    let tmp = x[jx as usize];
                    while i < min(n, j + k + 1) {
                        x[ix as usize] -= tmp * a[aj - j + i];
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
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let aj = j * lda;
            let mut i = j.saturating_sub(k);
            while i < j {
                tmp -= a[aj + k - j + i] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            if diag == 'n' || diag == 'N' {
                tmp /= a[aj + k];
            }
            x[jx as usize] = tmp;
            jx += incx;
            j += 1;
            if j > k {
                kx += incx;
            }
        }
    } else {
        kx += (n as isize - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let aj = j * lda;
            //FIXME Maybe we could do this in a diffrent way
            let mut i = min(n - 1, j + k);
            while i > j {
                tmp -= a[aj - j + i] * x[ix as usize];
                ix -= incx;
                i -= 1;
            }
            //FIXME
            if diag == 'n' || diag == 'N' {
                tmp /= a[aj];
            }
            x[jx as usize] = tmp;
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
    n: usize,
    ap: &[T],
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut kk = 0;
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                j += 1;
                let mut k = kk;
                kk += j;
                if !x[jx as usize].is_zero() {
                    let tmp = x[jx as usize];
                    let mut ix = kx;
                    while k < kk - 1 {
                        x[ix as usize] += tmp * ap[k];
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        x[jx as usize] *= ap[kk - 1];
                    }
                }
                jx += incx;
            }
        } else {
            let mut kk = n * (n + 1) / 2;
            kx += (n as isize - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                //FIXME figure out a better way
                if !x[jx as usize].is_zero() {
                    let tmp = x[jx as usize];
                    let mut ix = kx;
                    let mut k = kk - 1;
                    while k > kk - (n - j) {
                        x[ix as usize] += tmp * ap[k];
                        ix -= incx;
                        k -= 1;
                    }
                    if nounit {
                        x[jx as usize] *= ap[kk - n + j];
                    }
                }
                jx -= incx;
                kk -= n - j;
                //
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut kk = n * (n + 1) / 2;
        let mut jx = kx + (n as isize - 1) * incx;
        let mut j = n;
        while j >= 1 {
            kk -= j;
            j -= 1;
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let mut k = kk + j;
            if nounit {
                tmp *= ap[k];
            }
            while k > kk {
                ix -= incx;
                k -= 1;
                tmp += ap[k] * x[ix as usize];
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    } else {
        let mut kk = 0;
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let mut k = kk;
            kk += n - j;
            if nounit {
                tmp *= ap[k];
            }
            while k < kk {
                ix += incx;
                k += 1;
                tmp += ap[k] * x[ix as usize];
            }
            x[jx as usize] = tmp;
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
    n: usize,
    ap: &[T],
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut kk = n * (n + 1) / 2;
            let mut jx = kx + (n as isize - 1) * incx;
            let mut j = n;
            while j >= 1 {
                kk -= j;
                j -= 1;
                if !x[jx as usize].is_zero() {
                    let mut k = kk + j;
                    if nounit {
                        x[jx as usize] /= ap[k];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    while k > kk {
                        ix -= incx;
                        k -= 1;
                        x[ix as usize] -= tmp * ap[k];
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
                if !x[jx as usize].is_zero() {
                    if nounit {
                        x[jx as usize] /= ap[k];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    while k < kk {
                        k += 1;
                        ix += incx;
                        x[ix as usize] -= tmp * ap[k];
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
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let mut k = kk;
            while k < kk + j - 1 {
                tmp -= ap[k] * x[ix as usize];
                ix += incx;
                k += 1;
            }
            if nounit {
                // NOTE k == kk + j - 1
                tmp /= ap[k];
            }
            x[jx as usize] = tmp;
            jx += incx;
            kk += j;
        }
    } else {
        let mut kk = n * (n + 1) / 2;
        kx += (n as isize - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let mut k = kk - 1;
            kk -= n - j;
            while k > kk {
                tmp -= ap[k] * x[ix as usize];
                ix -= incx;
                k -= 1;
            }
            if nounit {
                // NOTE k == kk - n + j
                tmp /= ap[k];
            }
            x[jx as usize] = tmp;
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
    n: usize,
    a: &[T],
    lda: usize,
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = 0;
                    while i < j {
                        x[ix as usize] += tmp * a[aj + i];
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[aj + j];
                    }
                }
                jx += incx;
                j += 1;
            }
        } else {
            kx += (n as isize - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let aj = j * lda;
                    let mut ix = kx;
                    let mut i = n;
                    while i > j + 1 {
                        i -= 1;
                        x[ix as usize] += tmp * a[aj + i];
                        ix -= incx;
                    }
                    if nounit {
                        x[jx as usize] *= a[aj + j];
                    }
                }
                jx -= incx;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut jx = kx + (n as isize - 1) * incx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let aj = j * lda;
            let mut i = j;
            if nounit {
                tmp *= a[aj + j];
            }
            while i >= 1 {
                i -= 1;
                ix -= incx;
                tmp += a[aj + i] * x[ix as usize];
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let aj = j * lda;
            if nounit {
                tmp *= a[aj + j];
            }
            let mut i = j + 1;
            while i < n {
                ix += incx;
                tmp += a[aj + i] * x[ix as usize];
                i += 1;
            }
            x[jx as usize] = tmp;
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
    n: usize,
    a: &[T],
    lda: usize,
    x: &mut [T],
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
        kx = (-(n as isize) * incx) + incx
    };

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut jx = kx + (n as isize - 1) * incx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                if !x[jx as usize].is_zero() {
                    let aj = j * lda;
                    if nounit {
                        x[jx as usize] /= a[aj + j];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut i = j;
                    while i >= 1 {
                        i -= 1;
                        ix -= incx;
                        x[ix as usize] -= tmp * a[aj + i];
                    }
                }
                jx -= incx;
            }
        } else {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                if !x[jx as usize].is_zero() {
                    let aj = j * lda;
                    if nounit {
                        x[jx as usize] /= a[aj + j];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut i = j + 1;
                    while i < n {
                        ix += incx;
                        x[ix as usize] -= tmp * a[aj + i];
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
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                tmp -= a[aj + i] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            if nounit {
                tmp /= a[aj + i];
            }
            x[jx as usize] = tmp;
            jx += incx;
            j += 1;
        }
    } else {
        kx += (n as isize - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let aj = j * lda;
            let mut i = n;
            while i > j + 1 {
                i -= 1;
                tmp -= a[aj + i] * x[ix as usize];
                ix -= incx;
            }
            if nounit {
                tmp /= a[aj + j];
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    }
}
