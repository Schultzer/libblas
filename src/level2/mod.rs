// use std::arch::x86_64::*;
use num_traits::{Float, NumAssignOps};
use std::cmp::{max, min};

pub mod complex;

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
    a: &[T],
    lda: isize,
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
        kx = (-lenx * incx) + incx
    };
    if incy < 0 {
        ky = (-leny * incy) + incy
    };

    if !beta.is_one() {
        let mut iy = ky;
        let mut i = 0;
        while i < leny {
            y[iy as usize] = if beta.is_zero() {
                T::zero()
            } else {
                beta * y[iy as usize]
            };
            iy += incy;
            i += 1;
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
            let k = ku - j;
            let coor_aj = j * lda;
            let mut i = max(0, j - ku);
            while i <= min(m - 1, j + kl) {
                y[iy as usize] += tmp * a[(coor_aj + k + i) as usize];
                iy += incy;
                i += 1;
            }
            jx += incx;
            if j > ku - 1 {
                ky += incy
            };
            j += 1;
        }
    } else {
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let mut tmp = T::zero();
            let mut ix = kx;
            let k = ku - j;
            let coor_aj = j * lda;
            let mut i = max(0, j - ku);
            while i <= min(m - 1, j + kl) {
                tmp += a[(coor_aj + k + i) as usize] * x[ix as usize];
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
    m: isize,
    n: isize,
    alpha: T,
    a: &[T],
    lda: isize,
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
            let mut i = 0;
            while i < leny {
                y[i as usize] = T::zero();
                i += 1;
            }
        } else {
            let mut iy = ky;
            let mut i = 0;
            while i < leny {
                y[iy as usize] = if beta.is_zero() {
                    T::zero()
                } else {
                    y[iy as usize] * beta
                };
                iy += incy;
                i += 1;
            }
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
            let coor_aj = j * lda;
            let mut i = 0;
            while i < m {
                y[iy as usize] += tmp * a[(coor_aj + i) as usize];
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
            let coor_aj = j * lda;
            let mut i = 0;
            while i < m {
                tmp += a[(coor_aj + i) as usize] * x[ix as usize];
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
    m: isize,
    n: isize,
    alpha: T,
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    a: &mut [T],
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
        let mut tmp = y[jy as usize];
        if !tmp.is_zero() {
            tmp *= alpha;
            let mut ix = kx;
            let coords = j * lda;
            let mut i = 0;
            while i < m {
                a[(coords + i) as usize] += x[ix as usize] * tmp;
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
    a: &[T],
    lda: isize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
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
        let mut iy = ky;
        let mut i = 1;
        while i <= n {
            y[iy as usize] = if beta.is_zero() {
                T::zero()
            } else {
                beta * y[iy as usize]
            };
            iy += incy;
            i += 1;
        }
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
            let l = k - j;
            let coor_aj = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                y[iy as usize] += tmp * a[(coor_aj + l + i) as usize];
                tmp2 += a[(coor_aj + l + i) as usize] * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[(coor_aj + k) as usize] + alpha * tmp2;
            jx += incx;
            jy += incy;
            if j > k - 1 {
                kx += incx;
                ky += incy;
            }
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let coor_aj = j * lda;
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            y[jy as usize] += tmp * a[(coor_aj) as usize];
            let l = 1 - j;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i <= min(n - 1, j + k) {
                ix += incx;
                iy += incy;
                let index = ((coor_aj + l + i) - 1) as usize; // pretty hacky maybe there is a better way?
                y[iy as usize] += tmp * a[index];
                tmp2 += a[index] * x[ix as usize];
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
    n: isize,
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
    } else if n < 0 {
        info = 2;
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
        let mut iy = ky;
        let mut i = 0;
        while i < n {
            y[iy as usize] = if beta.is_zero() {
                T::zero()
            } else {
                beta * y[iy as usize]
            };
            iy += incy;
            i += 1;
        }
    }
    if alpha.is_zero() {
        return;
    }
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 1;
        while j <= n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            while k <= kk + j - 2 {
                y[iy as usize] += tmp * ap[k as usize];
                tmp2 += ap[k as usize] * x[ix as usize];
                ix += incx;
                iy += incy;
                k += 1;
            }
            y[jy as usize] += tmp * ap[(kk + j - 1) as usize] + alpha * tmp2;
            jx += incx;
            jy += incy;
            kk += j;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 1;
        while j <= n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            y[jy as usize] += tmp * ap[kk as usize];
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk + 1;
            while k <= kk + n - j {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * ap[k as usize];
                tmp2 += ap[k as usize] * x[ix as usize];
                k += 1;
            }
            y[jy as usize] += alpha * tmp2;
            jx += incx;
            jy += incy;
            kk += n - j + 1;
            j += 1;
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
    x: &[T],
    incx: isize,
    ap: &mut [T],
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
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
        let mut j = 1;
        while j <= n {
            if x[jx as usize] != T::zero() {
                let tmp = alpha * x[jx as usize];
                let mut ix = kx;
                let mut k = kk;
                while k < kk + j {
                    ap[k as usize] += x[ix as usize] * tmp;
                    ix += incx;
                    k += 1;
                }
            }
            jx += incx;
            kk += j;
            j += 1;
        }
    } else {
        let mut j = 1;
        while j <= n {
            let tmp = alpha * x[jx as usize];
            let mut ix = jx;
            let mut k = kk;
            while k <= kk + n - j {
                ap[k as usize] += x[ix as usize] * tmp;
                ix += incx;
                k += 1;
            }
            jx += incx;
            kk += n - j + 1;
            j += 1;
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
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    ap: &mut [T],
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
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
        let mut j = 1;
        while j <= n {
            let tmp = alpha * y[jy as usize];
            let tmp2 = alpha * x[jx as usize];
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            while k < kk + j {
                ap[k as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                ix += incx;
                iy += incy;
                k += 1;
            }
            jx += incx;
            jy += incy;
            kk += j;
            j += 1;
        }
    } else {
        let mut j = 1;
        while j <= n {
            let mut tmp = y[jy as usize];
            let mut tmp2 = x[jx as usize];
            if !tmp.is_zero() || tmp2.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = jx;
                let mut iy = jy;
                let mut k = kk;
                while k <= kk + n - j {
                    ap[k as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    ix += incx;
                    iy += incy;
                    k += 1;
                }
            }
            jx += incx;
            jy += incy;
            kk += n - j + 1;
            j += 1;
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
    a: &[T],
    lda: isize,
    x: &[T],
    incx: isize,
    beta: T,
    y: &mut [T],
    incy: isize,
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if n < 0 {
        info = 2;
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
        let mut iy = ky;
        let mut i = 1;
        while i <= n {
            y[iy as usize] = if beta.is_zero() {
                T::zero()
            } else {
                beta * y[iy as usize]
            };
            iy += incy;
            i += 1;
        }
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
            let coords = j * lda;
            let mut i = 0;
            while i < j {
                let index = (coords + i) as usize;
                y[iy as usize] += tmp * a[index];
                tmp2 += a[index] * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[(coords + j) as usize] + alpha * tmp2;
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2 = T::zero();
            let coords = j * lda;
            y[jy as usize] += tmp * a[(coords + j) as usize];
            let _l = 1 - j;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                let index = (coords + i) as usize;
                y[iy as usize] += tmp * a[index];
                tmp2 += a[index] * x[ix as usize];
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
    n: isize,
    alpha: T,
    x: &[T],
    incx: isize,
    a: &mut [T],
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
            let mut tmp = x[jx as usize];
            if !tmp.is_zero() {
                tmp *= alpha;
                let mut ix = kx;
                let coords = j * lda;
                let mut i = 0;
                while i <= j {
                    a[(coords + i) as usize] += x[ix as usize] * tmp;
                    ix += incx;
                    i += 1;
                }
            }
            jx += incx;
            j += 1;
        }
    } else {
        while j < n {
            let mut tmp = x[jx as usize];
            if !tmp.is_zero() {
                tmp *= alpha;
                let coords = j * lda;
                let mut ix = jx;
                let mut i = j;
                while i < n {
                    let delta = x[ix as usize] * tmp;
                    a[(coords + i) as usize] += delta;
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
    x: &[T],
    incx: isize,
    y: &[T],
    incy: isize,
    a: &mut [T],
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
            let mut tmp = y[jy as usize];
            let mut tmp2 = x[jx as usize];
            if !tmp2.is_zero() || !tmp.is_zero() {
                tmp *= alpha;
                tmp2 *= alpha;
                let mut ix = kx;
                let mut iy = ky;
                let coords = j * lda;
                let mut i = 0;
                while i <= j {
                    a[(coords + i) as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
                let coords = j * lda;
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j;
                while i < n {
                    a[(coords + i) as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
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
    a: &[T],
    lda: isize,
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let l = k - j;
                    let coords = j * lda;
                    let mut i = max(0, j - k);
                    while i < j {
                        x[ix as usize] += tmp * a[(coords + l + i) as usize];
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[(coords + k) as usize];
                    }
                }
                jx += incx;
                if j > k - 1 {
                    kx += incx;
                }
                j += 1;
            }
        } else {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n - 1;
            while j >= 0 {
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let l = 1 - j;
                    let coor_aj = j * lda;
                    let mut i = min(n - 1, j + k);
                    while i > j {
                        x[ix as usize] += tmp * a[((coor_aj + l + i) - 1) as usize]; // pretty hacky maybe there is a better way?
                        ix -= incx;
                        i -= 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[coor_aj as usize];
                    }
                }
                jx -= incx;
                if (n - 1) - j >= k {
                    kx -= incx;
                }
                j -= 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n - 1;
        while j >= 0 {
            let coor_aj = j * lda;
            let mut tmp = x[jx as usize];
            kx -= incx;
            let mut ix = kx;
            let l = k - j;
            if nounit {
                tmp *= a[(k + coor_aj) as usize];
            }
            let mut i = j - 1;
            while i >= max(0, j - k) {
                tmp += a[(l + i + coor_aj) as usize] * x[ix as usize];
                ix -= incx;
                i -= 1;
            }
            x[jx as usize] = tmp;
            jx -= incx;
            j -= 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let coor_aj = j * lda;
            let mut tmp = x[jx as usize];
            kx += incx;
            let mut ix = kx;
            let l = 1 - j;
            if nounit {
                tmp *= a[(coor_aj) as usize];
            }
            let mut i = j + 1;
            while i <= min(n - 1, j + k) {
                let index = (l + i + coor_aj - 1) as usize; // pretty hacky maybe there is a better way?
                tmp += a[index] * x[ix as usize];
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
    n: isize,
    k: isize,
    a: &[T],
    lda: isize,
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
            let mut j = n - 1;
            while j >= 0 {
                kx -= incx;
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let l = k - j;
                    let coords = j * lda;
                    if diag == 'n' || diag == 'N' {
                        x[jx as usize] /= a[(k + coords) as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut i = j - 1;
                    while i >= max(0, j - k) {
                        x[ix as usize] -= tmp * a[(coords + l + i) as usize];
                        ix -= incx;
                        i -= 1;
                    }
                }
                jx -= incx;
                j -= 1;
            }
        } else {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                kx += incx;
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let l = 1 - j;
                    let coords = j * lda;
                    if diag == 'n' || diag == 'N' {
                        x[jx as usize] /= a[(coords) as usize];
                    }
                    let mut i = j + 1;
                    let tmp = x[jx as usize];
                    while i <= min(n - 1, j + k) {
                        let index = ((coords + l + i) - 1) as usize; // Pretty hacky maybe there is a better way?
                        x[ix as usize] -= tmp * a[index];
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
            let l = k - j;
            let coords = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                tmp -= a[(coords + l + i) as usize] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            if diag == 'n' || diag == 'N' {
                tmp /= a[(coords + k) as usize];
            }
            x[jx as usize] = tmp;
            jx += incx;
            if j > k - 1 {
                kx += incx;
            }
            j += 1;
        }
    } else {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n - 1;
        while j >= 0 {
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let l = 1 - j;
            let coords = j * lda;
            let mut i = min(n - 1, j + k);
            while i > j {
                let index = ((l + i + coords) - 1) as usize; // Pretty hacky maybe there is a better way?
                tmp -= a[index] * x[ix as usize];
                ix -= incx;
                i -= 1;
            }
            if diag == 'n' || diag == 'N' {
                tmp /= a[(coords) as usize];
            }
            x[jx as usize] = tmp;
            jx -= incx;
            if (n - 1) - j >= k {
                kx -= incx;
            }
            j -= 1;
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
            let mut j = 1;
            while j <= n {
                if !x[jx as usize].is_zero() {
                    let tmp = x[jx as usize];
                    let mut ix = kx;
                    let mut k = kk;
                    while k <= kk + j - 2 {
                        x[ix as usize] += tmp * ap[k as usize];
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        x[jx as usize] *= ap[(kk + j - 1) as usize];
                    }
                }
                jx += incx;
                kk += j;
                j += 1;
            }
        } else {
            let mut kk = (n * (n + 1) / 2) - 1;
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                if !x[jx as usize].is_zero() {
                    let tmp = x[jx as usize];
                    let mut ix = kx;
                    let mut k = kk;
                    while k >= kk - (n - (j + 1)) {
                        x[ix as usize] += tmp * ap[k as usize];
                        ix -= incx;
                        k -= 1;
                    }
                    if nounit {
                        x[jx as usize] *= ap[(kk - n + j) as usize];
                    }
                }
                jx -= incx;
                kk -= n - j + 1;
                j -= 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut kk = ((n * (n + 1)) / 2) - 1;
        let mut jx = kx + (n - 1) * incx;
        let mut j = n;
        while j >= 1 {
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            if nounit {
                tmp *= ap[kk as usize];
            }
            let mut k = kk - 1;
            while k > kk - j {
                ix -= incx;
                tmp += ap[k as usize] * x[ix as usize];
                k -= 1;
            }
            x[jx as usize] = tmp;
            jx -= incx;
            kk -= j;
            j -= 1;
        }
    } else {
        let mut kk = 0;
        let mut jx = kx;
        let mut j = 1;
        while j <= n {
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            if nounit {
                tmp *= ap[kk as usize];
            }
            let mut k = kk + 1;
            while k <= kk + n - j {
                ix += incx;
                tmp += ap[k as usize] * x[ix as usize];
                k += 1;
            }
            x[jx as usize] = tmp;
            jx += incx;
            kk += n - j + 1;
            j += 1;
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
    } else if n < 0 {
        info = 4;
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
            let mut kk = ((n * (n + 1)) / 2) - 1;
            let mut jx = kx + (n - 1) * incx;
            let mut j = n;
            while j >= 1 {
                if !x[jx as usize].is_zero() {
                    if nounit {
                        x[jx as usize] /= ap[kk as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut k = kk - 1;
                    while k > kk - j {
                        ix -= incx;
                        x[ix as usize] -= tmp * ap[k as usize];
                        k -= 1;
                    }
                }
                jx -= incx;
                kk -= j;
                j -= 1;
            }
        } else {
            let mut kk = 0;
            let mut jx = kx;
            let mut j = 1;
            while j <= n {
                if !x[jx as usize].is_zero() {
                    if nounit {
                        x[jx as usize] /= ap[kk as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut k = kk + 1;
                    while k <= kk + n - j {
                        ix += incx;
                        x[ix as usize] -= tmp * ap[k as usize];
                        k += 1;
                    }
                }
                jx += incx;
                kk += n - j + 1;
                j += 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut kk = 0;
        let mut jx = kx;
        let mut j = 1;
        while j <= n {
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let mut k = kk;
            while k <= kk + j - 2 {
                tmp -= ap[k as usize] * x[ix as usize];
                ix += incx;
                k += 1;
            }
            if nounit {
                tmp /= ap[(kk + j - 1) as usize];
            }
            x[jx as usize] = tmp;
            jx += incx;
            kk += j;
            j += 1;
        }
    } else {
        let mut kk = (n * (n + 1) / 2) - 1;
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let mut k = kk;
            while k >= kk - (n - (j + 1)) {
                tmp -= ap[k as usize] * x[ix as usize];
                ix -= incx;
                k -= 1;
            }
            if nounit {
                tmp /= ap[(kk - n + j) as usize];
            }
            x[jx as usize] = tmp;
            jx -= incx;
            kk -= n - j + 1;
            j -= 1;
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
    a: &[T],
    lda: isize,
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let coords = j * lda;
                    let mut i = 0;
                    while i < j {
                        x[ix as usize] += tmp * a[(coords + i) as usize];
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[(coords + j) as usize];
                    }
                }
                jx += incx;
                j += 1;
            }
        } else {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n - 1;
            while j >= 0 {
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let coords = j * lda;
                    let mut ix = kx;
                    let mut i = n - 1;
                    while i > j {
                        x[ix as usize] += tmp * a[(coords + i) as usize];
                        ix -= incx;
                        i -= 1;
                    }
                    if nounit {
                        x[jx as usize] *= a[(coords + j) as usize];
                    }
                }
                jx -= incx;
                j -= 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut jx = kx + (n - 1) * incx;
        let mut j = n - 1;
        while j >= 0 {
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let coords = j * lda;
            let mut i = j - 1;
            if nounit {
                tmp *= a[(coords + j) as usize];
            }
            while i >= 0 {
                ix -= incx;
                tmp += a[(coords + i) as usize] * x[ix as usize];
                i -= 1;
            }
            x[jx as usize] = tmp;
            jx -= incx;
            j -= 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = x[jx as usize];
            let mut ix = jx;
            let coords = j * lda;
            if nounit {
                tmp *= a[(coords + j) as usize];
            }
            let mut i = j + 1;
            while i < n {
                ix += incx;
                tmp += a[(coords + i) as usize] * x[ix as usize];
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
    n: isize,
    a: &[T],
    lda: isize,
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
            let mut j = n - 1;
            while j >= 0 {
                if !x[jx as usize].is_zero() {
                    let coords = j * lda;
                    if nounit {
                        x[jx as usize] /= a[(coords + j) as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut i = j - 1;
                    while i >= 0 {
                        ix -= incx;
                        x[ix as usize] -= tmp * a[(coords + i) as usize];
                        i -= 1;
                    }
                }
                jx -= incx;
                j -= 1;
            }
        } else {
            let mut jx = kx;
            let mut j = 0;
            while j < n {
                if !x[jx as usize].is_zero() {
                    let coords = j * lda;
                    if nounit {
                        x[jx as usize] /= a[(coords + j) as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut ix = jx;
                    let mut i = j + 1;
                    while i < n {
                        ix += incx;
                        x[ix as usize] -= tmp * a[(coords + i) as usize];
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
            let coords = j * lda;
            let mut i = 0;
            while i < j {
                tmp -= a[(coords + i) as usize] * x[ix as usize];
                ix += incx;
                i += 1;
            }
            if nounit {
                tmp /= a[(coords + i) as usize];
            }
            x[jx as usize] = tmp;
            jx += incx;
            j += 1;
        }
    } else {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n - 1;
        while j >= 0 {
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            let coords = j * lda;
            let mut i = n - 1;
            while i > j {
                tmp -= a[(coords + i) as usize] * x[ix as usize];
                ix -= incx;
                i -= 1;
            }
            if nounit {
                tmp /= a[(coords + i) as usize];
            }
            x[jx as usize] = tmp;
            jx -= incx;
            j -= 1;
        }
    }
}
