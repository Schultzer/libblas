// use std::arch::x86_64::*;
use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::{max, min};

fn multiply<T: Float + NumAssignOps>(
    left: &mut [Complex<T>],
    right: Complex<T>,
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

fn zero<T: Float + NumAssignOps>(
    left: &mut [Complex<T>],
    len: isize,
    mut index: isize,
    inc: isize,
) {
    let mut i = 0;
    while i < len {
        left[index as usize] = Complex::zero();
        index += inc;
        i += 1;
    }
}

/// GBMV  performs one of the matrix-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or y := alpha*A**H*x + beta*y,
/// where alpha and beta are scalars, x and y are vectors and A is an m by n band matrix, with kl sub-diagonals and ku super-diagonals.
/// This is [CGBMV](http://www.netlib.org/lapack/explore-html/d0/d75/cgbmv_8f.html) and [ZGBMV](http://www.netlib.org/lapack/explore-html/d9/d46/zgbmv_8f.html) comined in one function
#[inline]
pub fn gbmv<T: Float + NumAssignOps>(
    trans: char,
    m: usize,
    n: usize,
    kl: usize,
    ku: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
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

    let noconj = trans == 't' || trans == 'T';

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
            let tmp = x[jx as usize] * alpha;
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
                ky += incy;
            };
            j += 1;
            //
        }
    } else {
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let mut tmp: Complex<T> = Complex::zero();
            let mut ix = kx;
            // FIXME
            let k = ku as isize - j as isize; // apparently we want this to be negative sometimes.
            let aj = j * lda;
            let mut i = j.saturating_sub(ku); // MAX(0, J-KU) using saturating_sub it will always be 0
            while i < min(m, j + kl + 1) {
                if noconj {
                    tmp += a[(aj as isize + k + i as isize) as usize] * x[ix as usize];
                } else {
                    tmp += a[(aj as isize + k + i as isize) as usize].conj() * x[ix as usize];
                }
                ix += incx;
                i += 1;
            }
            y[jy as usize] += alpha * tmp;
            jy += incy;
            if j > ku - 1 {
                kx += incx
            };
            j += 1;
            //
        }
    }
}

/// GEMV  performs one of the matrix-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or  y := alpha*A**H*x + beta*y,
/// where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
/// This is [CGEMV](http://www.netlib.org/lapack/explore-html/d4/d8a/cgemv_8f.html) and [ZGEMV](http://www.netlib.org/lapack/explore-html/db/d40/zgemv_8f.html) comined in one function
#[inline]
pub fn gemv<T: Float + NumAssignOps>(
    trans: char,
    m: usize,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
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

    if m == 0 || n == 0 || (alpha.is_zero() && beta.is_zero()) {
        return;
    }

    let noconj = trans == 't' || trans == 'T';
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
            let mut tmp: Complex<T> = Complex::zero();
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                if noconj {
                    tmp += a[aj + i] * x[ix as usize];
                } else {
                    tmp += a[aj + i].conj() * x[ix as usize];
                }

                ix += incx;
                i += 1;
            }
            y[jy as usize] += alpha * tmp;
            jy += incy;
            j += 1;
        }
    }
}

/// GERC  performs the rank 1 operation
/// A := alpha*x*y**H + A,
/// where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix.
/// This is [CGERC](http://www.netlib.org/lapack/explore-html/dd/d84/cgerc_8f.html) and [ZGERC](http://www.netlib.org/lapack/explore-html/d3/dad/zgerc_8f.html) comined in one function
#[inline]
pub fn gerc<T: Float + NumAssignOps>(
    m: usize,
    n: usize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
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
        panic!("gerc {}", info);
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
        let tmp = y[jy as usize];
        if !tmp.is_zero() {
            let tmp = alpha * tmp.conj();
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

/// GERU  performs the rank 1 operation
///  A := alpha*x*y**T + A,
/// where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix.
/// This is [CGERU](http://www.netlib.org/lapack/explore-html/db/d5f/cgeru_8f.html) and [ZGERU](http://www.netlib.org/lapack/explore-html/d7/d12/zgeru_8f.html) comined in one function
#[inline]
pub fn geru<T: Float + NumAssignOps>(
    m: usize,
    n: usize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
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
        panic!("geru {}", info);
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
        let tmp = y[jy as usize];
        if !tmp.is_zero() {
            let tmp = alpha * tmp;
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

/// HBMV  performs the matrix-vector  operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n hermitian band matrix, with k super-diagonals.
/// This is [CHBMV](http://www.netlib.org/lapack/explore-html/db/dc2/chbmv_8f.html) and [ZHBMV](http://www.netlib.org/lapack/explore-html/d3/d1a/zhbmv_8f.html) comined in one function
#[inline]
pub fn hbmv<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    k: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
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
        panic!("hbmv {}", info);
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
    if uplo == 'u' || uplo == 'U' {
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = j.saturating_sub(k); // MAX(0, j - k)
            while i < j {
                y[iy as usize] += tmp * a[aj + k - j + i];
                tmp2 += a[aj + k - j + i].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[iy as usize] += tmp * a[aj + k].re + (alpha * tmp2);
            jx += incx;
            jy += incy;
            j += 1;
            if j > k {
                kx += incx;
                ky += incy;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            let aj = j * lda;
            y[jy as usize] += tmp * a[aj].re;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[aj - j + i];
                tmp2 += a[aj - j + i].conj() * x[ix as usize];
                i += 1;
            }
            y[jy as usize] += alpha * tmp2;
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// HEMV  performs the matrix-vector  operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n hermitian matrix.
/// This is [CHEMV](http://www.netlib.org/lapack/explore-html/d7/d51/chemv_8f.html) and [ZHEMV](http://www.netlib.org/lapack/explore-html/d0/ddd/zhemv_8f.html) comined in one function
#[inline]
pub fn hemv<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
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
        panic!("hemv {}", info);
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
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                y[iy as usize] += tmp * a[aj + i];
                tmp2 += a[aj + i].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[aj + j].re + (alpha * tmp2);
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        while j < n {
            let aj = j * lda;
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            y[jy as usize] += tmp * a[aj + j].re;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[aj + i];
                tmp2 += a[aj + i].conj() * x[ix as usize];
                i += 1;
            }
            y[jy as usize] += alpha * tmp2;
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// HER   performs the hermitian rank 1 operation
///  A := alpha*x*x**H + A,
///  where alpha is a real scalar, x is an n element vector and A is an n by n hermitian matrix.
/// This is [CHER](http://www.netlib.org/lapack/explore-html/d3/d6d/cher_8f.html) and [ZHER](http://www.netlib.org/lapack/explore-html/de/d0e/zher_8f.html) comined in one function
#[inline]
pub fn her<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: T,
    x: &[Complex<T>],
    incx: isize,
    a: &mut [Complex<T>],
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
        panic!("her {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    if incx < 0 {
        kx = (-(n as isize) * incx) + incx
    };
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let tmp = x[jx as usize];
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let mut ix = kx;
                let mut i = 0;
                while i < j {
                    a[aj + i] += tmp * x[ix as usize];
                    ix += incx;
                    i += 1;
                }
                let tmp = tmp * x[jx as usize];
                a[aj + j] = Complex {
                    re: a[aj + j].re + tmp.re,
                    im: T::zero(),
                };
            } else {
                a[aj + j].im = T::zero();
            }
            jx += incx;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let tmp = x[jx as usize];
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let tmp2 = tmp * x[jx as usize];
                a[aj + j] = Complex {
                    re: a[aj + j].re + tmp2.re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    a[aj + i] += tmp * x[ix as usize];
                    i += 1;
                }
            } else {
                a[aj + j].im = T::zero();
            }
            jx += incx;
            j += 1;
        }
    }
}

/// HER2  performs the hermitian rank 2 operation
/// A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
/// where alpha is a scalar, x and y are n element vectors and A is an n by n hermitian matrix.
/// This is [CHER2](http://www.netlib.org/lapack/explore-html/db/d87/cher2_8f.html) and [ZHER2](http://www.netlib.org/lapack/explore-html/da/d8a/zher2_8f.html) comined in one function
#[inline]
pub fn her2<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
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
        panic!("her2 {}", info);
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
            let aj = j * lda;
            let y1 = y[jy as usize];
            let x1 = x[jx as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();

                let mut ix = kx;
                let mut iy = ky;
                let mut i = 0;
                while i < j {
                    a[aj + i] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    ix += incx;
                    iy += incy;
                    i += 1;
                }
                // FIXME
                // maybe only calculate the real.
                a[aj + j] = Complex {
                    re: a[aj + j].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
            } else {
                a[aj + j].im = T::zero()
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let aj = j * lda;
            let y1 = y[jy as usize];
            let x1 = x[jx as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                a[aj + j] = Complex {
                    re: a[aj + j].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    iy += incy;
                    a[aj + i] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    i += 1;
                }
            } else {
                a[aj + j].im = T::zero()
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    }
}

/// HPMV  performs the matrix-vector operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n hermitian matrix, supplied in packed form.
/// This is [CHPMV](http://www.netlib.org/lapack/explore-html/d2/d06/chpmv_8f.html) and [ZHPMV](http://www.netlib.org/lapack/explore-html/d0/d60/zhpmv_8f.html) comined in one function
#[inline]
pub fn hpmv<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: Complex<T>,
    ap: &[Complex<T>],
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
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
        panic!("hpmv {}", info);
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
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j;
            while k < kk - 1 { // KK + J - 2
                y[iy as usize] += tmp * ap[k];
                tmp2 += ap[k].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                k += 1;
            }
            y[jy as usize] += tmp * ap[kk - 1].re + alpha * tmp2;
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
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk;
            kk += n - j;
            y[jy as usize] += tmp * ap[k].re;
            while k < kk {
                k += 1;
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * ap[k];
                tmp2 += ap[k].conj() * x[ix as usize];
            }
            y[jy as usize] += alpha * tmp2;
            jx += incx;
            jy += incy;
            kk += 1;
        }
    }
}

/// HPR    performs the hermitian rank 1 operation
/// A := alpha*x*x**H + A,
/// where alpha is a real scalar, x is an n element vector and A is an n by n hermitian matrix, supplied in packed form.
/// This is [CHPR](http://www.netlib.org/lapack/explore-html/db/dcd/chpr_8f.html) and [ZHPR](http://www.netlib.org/lapack/explore-html/de/de1/zhpr_8f.html) comined in one function
#[inline]
pub fn hpr<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: T,
    x: &[Complex<T>],
    incx: isize,
    ap: &mut [Complex<T>],
) {
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if incx == 0 {
        info = 6;
    }
    if info != 0 {
        panic!("hpr {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }
    let mut kx = 0;
    if incx < 0 {
        kx = (-(n as isize) * incx) + incx
    };
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let x1 = x[jx as usize];
            let mut k = kk;
            kk += j;
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                let mut ix = kx;
                while k < kk {
                    ap[k] += tmp * x[ix as usize];
                    ix += incx;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                ap[kk] = Complex {
                    re: ap[kk].re + (x1 * tmp).re,
                    im: T::zero(),
                };
            } else {
                ap[kk].im = T::zero();
            }
            jx += incx;
            kk += 1;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            j += 1;
            let x1 = x[jx as usize];
            let mut k = kk;
            kk += n - j;
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                ap[k] = Complex {
                    re: ap[k].re + (x1 * tmp).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                while k < kk {
                    k += 1;
                    ix += incx;
                    ap[k] += tmp * x[ix as usize];
                }
            } else {
                ap[k].im = T::zero();
            }
            jx += incx;
            kk += 1;
        }
    }
}

/// HPR2  performs the hermitian rank 2 operation
/// A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
/// where alpha is a scalar, x and y are n element vectors and A is an n by n hermitian matrix, supplied in packed form.
/// This is [CHPR2](http://www.netlib.org/lapack/explore-html/d6/d44/chpr2_8f.html) and [ZHPR2](http://www.netlib.org/lapack/explore-html/d5/d52/zhpr2_8f.html) comined in one function
#[inline]
pub fn hpr2<T: Float + NumAssignOps>(
    uplo: char,
    n: usize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    ap: &mut [Complex<T>],
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
        panic!("hpr2 {}", info);
    }

    let alpha_is_zero = alpha.is_zero();

    if n == 0 || alpha_is_zero {
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
            let x1 = x[jx as usize];
            let y1 = y[jy as usize];
            let mut k = kk;
            kk += j;
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                let mut ix = kx;
                let mut iy = ky;
                while k < kk {
                    ap[k] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    ix += incx;
                    iy += incy;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                ap[kk] = Complex {
                    re: ap[kk].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
            } else {
                ap[kk].im = T::zero();
            }
            jx += incx;
            jy += incy;
            kk += 1;
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            j += 1;
            let x1 = x[jx as usize];
            let y1 = y[jy as usize];
            let mut k = kk;
            kk += n - j;
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                ap[k] = Complex {
                    re: ap[k].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut iy = jy;
                // FIXME Maybe we could move kk += n - j + 1 up here
                while k < kk {
                    k += 1;
                    ix += incx;
                    iy += incy;
                    ap[k] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                }
            } else {
                ap[k].im = T::zero();
            }
            jx += incx;
            jy += incy;
            kk += 1;
        }
    }
}

/// TBMV  performs one of the matrix-vector operations
/// x := A*x,   or   x := A**T*x,   or   x := A**H*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular band matrix, with ( k + 1 ) diagonals.
/// This is [CTBMV](http://www.netlib.org/lapack/explore-html/d3/dcd/ctbmv_8f.html) and [ZTBMV](http://www.netlib.org/lapack/explore-html/d3/d39/ztbmv_8f.html) comined in one function
#[inline]
pub fn tbmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    k: usize,
    a: &[Complex<T>],
    lda: usize,
    x: &mut [Complex<T>],
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
    let noconj = trans == 't' || trans == 'T';
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
                        x[ix as usize] += tmp * a[aj + i + k - j];
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        x[jx as usize] = tmp * a[aj + k]
                    }
                }
                jx += incx;
                j += 1;
                if j > k {
                    kx += incx;
                }
            }
        } else {
            // FIXME
            kx += (n as isize - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let aj = j * lda;
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let mut i = min(n - 1, j + k); // NOTE
                    while i >= j + 1 {
                        x[ix as usize] += tmp * a[aj - j + i];
                        ix -= incx;
                        i -= 1;
                    }
                    if nounit {
                        x[jx as usize] = tmp * a[aj];
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
            let mut tmp = x[jx as usize];
            kx -= incx;
            let mut ix = kx;
            let aj = j * lda;
            if nounit {
                if noconj {
                    tmp *= a[aj + k];
                } else {
                    tmp *= a[aj + k].conj();
                }
            }
            let mut i = j;
            // FIXME figure out if there is a better way
            while i >= max(1, j.saturating_sub(k) + 1) {
                i -= 1;
                if noconj {
                    tmp += a[aj + i + k - j] * x[ix as usize];
                } else {
                    tmp += a[aj + i + k - j].conj() * x[ix as usize];
                }
                ix -= incx;
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = x[jx as usize];
            kx += incx;
            let aj = j * lda;
            if nounit {
                if noconj {
                    tmp *= a[aj];
                } else {
                    tmp *= a[aj].conj();
                }
            }
            let mut ix = kx;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                if noconj {
                    tmp += a[aj - j + i] * x[ix as usize];
                } else {
                    tmp += a[aj - j + i].conj() * x[ix as usize];
                }

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
/// A*x = b,   or   A**T*x = b,   or   A**H*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular band matrix, with ( k + 1 )  diagonals.
/// No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.
/// This is [CTBSV](http://www.netlib.org/lapack/explore-html/d9/d5f/ctbsv_8f.html) and [ZTBSV](http://www.netlib.org/lapack/explore-html/d4/d5a/ztbsv_8f.html) comined in one function
#[inline]
pub fn tbsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    k: usize,
    a: &[Complex<T>],
    lda: usize,
    x: &mut [Complex<T>],
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

    let noconj = trans == 't' || trans == 'T';
    let nounit = diag == 'n' || diag == 'N';
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
                    if nounit {
                        x[jx as usize] /= a[aj + k];
                    }
                    let tmp = x[jx as usize];
                    let mut i = j;
                    while i > j.saturating_sub(k) {
                        i -= 1;
                        x[ix as usize] -= tmp * a[aj +  k - j + i];
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
                    if nounit {
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
                if noconj {
                    tmp -= a[aj + k - j + i] * x[ix as usize];
                } else {
                    tmp -= a[aj + k - j + i].conj() * x[ix as usize];
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[aj + k];
                } else {
                    tmp /= a[aj + k].conj();
                };
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
            let aj = j * lda;
            let mut tmp = x[jx as usize];
            let mut ix = kx;
            //FIXME Maybe we could do this in a diffrent way
            let mut i = min(n - 1, j + k);
            while i > j {
                if noconj {
                    tmp -= a[aj - j + i] * x[ix as usize];
                } else {
                    tmp -= a[aj - j + i].conj() * x[ix as usize];
                };
                i -= 1;
                ix -= incx;
            }
            //
            if nounit {
                if noconj {
                    tmp /= a[aj];
                } else {
                    tmp /= a[aj].conj();
                };
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
///  x := A*x,   or   x := A**T*x,   or   x := A**H*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// This is [CTPMV](http://www.netlib.org/lapack/explore-html/d4/dbb/ctpmv_8f.html) and [ZTPMV](http://www.netlib.org/lapack/explore-html/d2/d9e/ztpmv_8f.html) comined in one function
#[inline]
pub fn tpmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    ap: &[Complex<T>],
    x: &mut [Complex<T>],
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

    let noconj = trans == 't' || trans == 'T';
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
                if !x[jx as usize].is_zero() {
                    let tmp = x[jx as usize];
                    let mut ix = kx;
                    let mut k = kk;
                    while k < kk + j - 1 {
                        x[ix as usize] += tmp * ap[k];
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        x[jx as usize] *= ap[kk + j - 1];
                        // x[jx as usize] = tmp * ap[kk + j - 1];
                    }
                }
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
                if noconj {
                    tmp *= ap[k];
                } else {
                    tmp *= ap[k].conj();
                };
            }
            while k > kk {
                ix -= incx;
                k -= 1;
                if noconj {
                    tmp += ap[k] * x[ix as usize];
                } else {
                    tmp += ap[k].conj() * x[ix as usize];
                };
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
                if noconj {
                    tmp *= ap[k];
                } else {
                    tmp *= ap[k].conj();
                };
            }
            while k < kk {
                ix += incx;
                k += 1;
                if noconj {
                    tmp += ap[k] * x[ix as usize];
                } else {
                    tmp += ap[k].conj() * x[ix as usize];
                };
            }
            x[jx as usize] = tmp;
            jx += incx;
            kk += 1;
        }
    }
}

/// TPSV  solves one of the systems of equations
/// A*x = b,   or   A**T*x = b,   or   A**H*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.
/// This is [CTPSV](http://www.netlib.org/lapack/explore-html/d8/d56/ctpsv_8f.html) and [ZTPSV](http://www.netlib.org/lapack/explore-html/da/d57/ztpsv_8f.html) comined in one function
#[inline]
pub fn tpsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    ap: &[Complex<T>],
    x: &mut [Complex<T>],
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
    let noconj = trans == 't' || trans == 'T';
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
            kk += j;
            while k < kk - 1 {
                if noconj {
                    tmp -= ap[k] * x[ix as usize];
                } else {
                    tmp -= ap[k].conj() * x[ix as usize];
                };
                ix += incx;
                k += 1;
            }
            if nounit {
                // NOTE k == kk + j - 1
                if noconj {
                    tmp /= ap[k];
                } else {
                    tmp /= ap[k].conj();
                };
            }
            x[jx as usize] = tmp;
            jx += incx;
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
                if noconj {
                    tmp -= ap[k] * x[ix as usize];
                } else {
                    tmp -= ap[k].conj() * x[ix as usize];
                };

                ix -= incx;
                k -= 1;
            }
            if nounit {
                // NOTE k == kk - n + j
                if noconj {
                    tmp /= ap[k];
                } else {
                    tmp /= ap[k].conj();
                };
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    }
}

/// TRMV  performs one of the matrix-vector operations
/// x := A*x,   or   x := A**T*x,   or   x := A**H*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix.
/// This is [CTRMV](http://www.netlib.org/lapack/explore-html/df/d78/ctrmv_8f.html) and [ZTRMV](http://www.netlib.org/lapack/explore-html/d0/dd1/ztrmv_8f.html) comined in one function
#[inline]
pub fn trmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    a: &[Complex<T>],
    lda: usize,
    x: &mut [Complex<T>],
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

    let noconj = trans == 't' || trans == 'T';
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
                        x[jx as usize] = tmp * a[aj + j];
                    }
                }
                j += 1;
                jx += incx;
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
                        x[jx as usize] = tmp * a[aj + j];
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
            if noconj {
                if nounit {
                    tmp *= a[aj + j];
                }
                while i >= 1 {
                    i -= 1;
                    ix -= incx;
                    tmp += a[aj + i] * x[ix as usize];
                }
            } else {
                if nounit {
                    tmp *= a[aj + j].conj();
                }
                while i >= 1 {
                    i -= 1;
                    ix -= incx;
                    tmp += a[aj + i].conj() * x[ix as usize];
                }
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
            let mut i = j + 1;
            if noconj {
                if nounit {
                    tmp *= a[aj + j];
                }
                while i < n {
                    ix += incx;
                    tmp += a[aj + i] * x[ix as usize];
                    i += 1;
                }
            } else {
                if nounit {
                    tmp *= a[aj + j].conj();
                }
                while i < n {
                    ix += incx;
                    tmp += a[aj + i].conj() * x[ix as usize];
                    i += 1;
                }
            }
            x[jx as usize] = tmp;
            jx += incx;
            j += 1;
        }
    }
}

/// TRSV  solves one of the systems of equations
/// A*x = b,   or   A**T*x = b,   or   A**H*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular matrix.
/// No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.
/// This is [CTRSV](http://www.netlib.org/lapack/explore-html/d4/dc8/ctrsv_8f.html) and [ZTRSV](http://www.netlib.org/lapack/explore-html/d1/d2f/ztrsv_8f.html) comined in one function
#[inline]
pub fn trsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: usize,
    a: &[Complex<T>],
    lda: usize,
    x: &mut [Complex<T>],
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

    let noconj = trans == 't' || trans == 'T';
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let aj = j * lda;
                    if nounit {
                        x[jx as usize] = tmp / a[aj + j];
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let aj = j * lda;
                    if nounit {
                        x[jx as usize] = tmp / a[aj + j];
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
            let mut ix = kx;
            let mut tmp = x[jx as usize];
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                if noconj {
                    tmp -= a[aj + i] * x[ix as usize];
                } else {
                    tmp -= a[aj + i].conj() * x[ix as usize];
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[aj + j];
                } else {
                    tmp /= a[aj + j].conj();
                };
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
            let mut ix = kx;
            let mut tmp = x[jx as usize];
            let aj = j * lda;
            let mut i = n;
            while i > j + 1 {
                i -= 1;
                if noconj {
                    tmp -= a[aj + i] * x[ix as usize];
                } else {
                    tmp -= a[aj + i].conj() * x[ix as usize];
                };
                ix -= incx;
            }
            if nounit {
                if noconj {
                    tmp /= a[aj + j];
                } else {
                    tmp /= a[aj + j].conj();
                };
            }
            x[jx as usize] = tmp;
            jx -= incx;
        }
    }
}
