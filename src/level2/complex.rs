// use std::arch::x86_64::*;
use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::{max, min};

/// GBMV  performs one of the matrix-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or y := alpha*A**H*x + beta*y,
/// where alpha and beta are scalars, x and y are vectors and A is an m by n band matrix, with kl sub-diagonals and ku super-diagonals.
/// This is [CGBMV](http://www.netlib.org/lapack/explore-html/d0/d75/cgbmv_8f.html) and [ZGBMV](http://www.netlib.org/lapack/explore-html/d9/d46/zgbmv_8f.html) comined in one function
#[inline]
pub fn gbmv<T: Float + NumAssignOps>(
    trans: char,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: isize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
    incy: isize,
) {
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();
    let alpha_is_zero = alpha.is_one();

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

    if m == 0 || n == 0 || (alpha_is_zero && beta_is_zero) {
        return;
    }

    let lenx = if trans == 'n' || trans == 'N' { n } else { m };
    let leny = if trans == 'n' || trans == 'N' { m } else { n };

    let noconj = trans == 't' || trans == 'T';

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-lenx * incx) + incx
    };
    if incy < 0 {
        ky = (-leny * incy) + incy
    };

    if !beta_is_one {
        let mut iy = ky;
        let mut i = 1;
        while i <= leny {
            if beta_is_zero {
                y[iy as usize] = Complex::zero();
            } else {
                y[iy as usize] *= beta;
            }
            iy += incy;
            i += 1;
        }
    }

    if alpha_is_zero {
        return;
    }
    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = x[jx as usize] * alpha;
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
                ky += incy;
            };
            j += 1;
        }
    } else {
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            let mut tmp: Complex<T> = Complex::zero();
            let mut ix = kx;
            let k = ku - j;
            let coor_aj = j * lda;
            let istart = max(0, j - ku);
            let iend = min(m - 1, j + kl);
            let mut i = istart;
            while i <= iend {
                if noconj {
                    tmp += a[(coor_aj + k + i) as usize] * x[ix as usize];
                } else {
                    tmp += a[(coor_aj + k + i) as usize].conj() * x[ix as usize];
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
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: isize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
    incy: isize,
) {
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();
    let alpha_is_zero = alpha.is_zero();
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

    if m == 0 || n == 0 || (alpha_is_zero && beta_is_zero) {
        return;
    }

    let noconj = trans == 't' || trans == 'T';

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

    if !beta_is_one {
        let mut iy = ky;
        let mut i = 0;
        while i < leny {
            if beta_is_zero {
                y[iy as usize] = Complex::zero();
            } else {
                y[iy as usize] *= beta;
            }
            iy += incy;
            i += 1;
        }
    }
    if alpha_is_zero {
        return;
    }
    if trans == 'n' || trans == 'N' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut iy = ky;
            let coords = j * lda;
            let mut i = 0;
            while i < m {
                y[iy as usize] += tmp * a[(coords + i) as usize];
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
            let coords = j * lda;
            let mut i = 0;
            while i < m {
                if noconj {
                    tmp += a[(coords + i) as usize] * x[ix as usize];
                } else {
                    tmp += a[(coords + i) as usize].conj() * x[ix as usize];
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
    m: isize,
    n: isize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
    lda: isize,
) {
    let alpha_is_zero = alpha.is_zero();
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
        panic!("gerc {}", info);
    }

    if m == 0 || n == 0 || alpha_is_zero {
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
        let tmp = y[jy as usize];
        if !tmp.is_zero() {
            let tmp = alpha * tmp.conj();
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

/// GERU  performs the rank 1 operation
///  A := alpha*x*y**T + A,
/// where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix.
/// This is [CGERU](http://www.netlib.org/lapack/explore-html/db/d5f/cgeru_8f.html) and [ZGERU](http://www.netlib.org/lapack/explore-html/d7/d12/zgeru_8f.html) comined in one function
#[inline]
pub fn geru<T: Float + NumAssignOps>(
    m: isize,
    n: isize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
    lda: isize,
) {
    let alpha_is_zero = alpha.is_zero();
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
        panic!("geru {}", info);
    }

    if m == 0 || n == 0 || alpha_is_zero {
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
        let tmp = y[jy as usize];
        if !tmp.is_zero() {
            let tmp = alpha * tmp;
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

/// HBMV  performs the matrix-vector  operation
/// y := alpha*A*x + beta*y,
/// where alpha and beta are scalars, x and y are n element vectors and A is an n by n hermitian band matrix, with k super-diagonals.
/// This is [CHBMV](http://www.netlib.org/lapack/explore-html/db/dc2/chbmv_8f.html) and [ZHBMV](http://www.netlib.org/lapack/explore-html/d3/d1a/zhbmv_8f.html) comined in one function
#[inline]
pub fn hbmv<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    k: isize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: isize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
    incy: isize,
) {
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();
    let alpha_is_zero = alpha.is_zero();
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
        panic!("hbmv {}", info);
    }

    if n == 0 || (alpha_is_zero && beta_is_one) {
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

    if !beta_is_one {
        let mut iy = ky;
        let mut i = 0;
        while i < n {
            if beta_is_zero {
                y[iy as usize] = Complex::zero();
            } else {
                y[iy as usize] *= beta;
            }
            iy += incy;
            i += 1;
        }
    }
    if alpha_is_zero {
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
            let l = k - j;
            let coords = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                y[iy as usize] += tmp * a[(coords + l + i) as usize];
                tmp2 += a[(coords + l + i) as usize].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[iy as usize] += tmp * a[(coords + k) as usize].re + (alpha * tmp2);
            jx += incx;
            jy += incy;
            if j > k - 1 {
                kx += incx;
                ky += incy;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            let coords = j * lda;
            y[jy as usize] += tmp * a[coords as usize].re;
            let l = 1 - j;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i <= min(n - 1, j + k) {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[((coords + l + i) - 1) as usize];
                tmp2 += a[((coords + l + i) - 1) as usize].conj() * x[ix as usize];
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
    n: isize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: isize,
    x: &[Complex<T>],
    incx: isize,
    beta: Complex<T>,
    y: &mut [Complex<T>],
    incy: isize,
) {
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();
    let alpha_is_zero = alpha.is_zero();
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
        panic!("hemv {}", info);
    }

    if n == 0 || (alpha_is_zero && beta_is_one) {
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

    if !beta_is_one {
        let mut iy = ky;
        let mut i = 0;
        while i < n {
            if beta_is_zero {
                y[iy as usize] = Complex::zero();
            } else {
                y[iy as usize] *= beta;
            }
            iy += incy;
            i += 1;
        }
    }
    if alpha_is_zero {
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
            let coords = j * lda;
            let mut i = 0;
            while i < j {
                y[iy as usize] += tmp * a[(coords + i) as usize];
                tmp2 += a[(coords + i) as usize].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                i += 1;
            }
            y[jy as usize] += tmp * a[(coords + j) as usize].re + (alpha * tmp2);
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        while j < n {
            let coords = j * lda;
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            y[jy as usize] += tmp * a[(coords + j) as usize].re;
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * a[(coords + i) as usize];
                tmp2 += a[(coords + i) as usize].conj() * x[ix as usize];
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
    n: isize,
    alpha: T,
    x: &[Complex<T>],
    incx: isize,
    a: &mut [Complex<T>],
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
        panic!("her {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }

    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let coords = j * lda;
            let tmp = x[jx as usize];
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let mut ix = kx;
                let mut i = 0;
                while i < j {
                    a[(coords + i) as usize] += tmp * x[ix as usize];
                    ix += incx;
                    i += 1;
                }
                let tmp = tmp * x[jx as usize];
                a[(coords + j) as usize] = Complex {
                    re: a[(coords + j) as usize].re + tmp.re,
                    im: T::zero(),
                };
            } else {
                a[(coords + j) as usize].im = T::zero();
            }
            jx += incx;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let coords = j * lda;
            let tmp = x[jx as usize];
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let tmp2 = tmp * x[jx as usize];
                a[(coords + j) as usize] = Complex {
                    re: a[(coords + j) as usize].re + tmp2.re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    a[(coords + i) as usize] += tmp * x[ix as usize];
                    i += 1;
                }
            } else {
                a[(coords + j) as usize].im = T::zero();
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
    n: isize,
    alpha: Complex<T>,
    x: &[Complex<T>],
    incx: isize,
    y: &[Complex<T>],
    incy: isize,
    a: &mut [Complex<T>],
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
        panic!("her2 {}", info);
    }

    let alpha_is_zero = alpha.is_zero();

    if n == 0 || alpha_is_zero {
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
            let coords = j * lda;
            let y1 = y[jy as usize];
            let x1 = x[jx as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();

                let mut ix = kx;
                let mut iy = ky;
                let mut i = 0;
                while i < j {
                    a[(coords + i) as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    ix += incx;
                    iy += incy;
                    i += 1;
                }
                // FIXME
                // maybe only calculate the real.
                a[(coords + j) as usize] = Complex {
                    re: a[(coords + j) as usize].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
            } else {
                a[(coords + j) as usize].im = T::zero()
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let coords = j * lda;
            let y1 = y[jy as usize];
            let x1 = x[jx as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                a[(coords + j) as usize] = Complex {
                    re: a[(coords + j) as usize].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    iy += incy;
                    a[(coords + i) as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    i += 1;
                }
            } else {
                a[(coords + j) as usize].im = T::zero()
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
    n: isize,
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
    } else if n < 0 {
        info = 2;
    } else if incx == 0 {
        info = 6;
    } else if incy == 0 {
        info = 9;
    }
    if info != 0 {
        panic!("hpmv {}", info);
    }

    let beta_is_one = beta.is_one();
    let beta_is_zero = beta.is_zero();
    let alpha_is_zero = alpha.is_zero();

    if n == 0 || (alpha_is_zero && beta_is_one) {
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

    if !beta_is_one {
        let mut iy = ky;
        let mut i = 0;
        while i < n {
            if beta_is_zero {
                y[iy as usize] = Complex::zero();
            } else {
                y[iy as usize] *= beta;
            }
            iy += incy;
            i += 1;
        }
    }
    if alpha_is_zero {
        return;
    }
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 1;
        while j <= n {
            let tmp = alpha * x[jx as usize];
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            while k <= kk + j - 2 {
                y[iy as usize] += tmp * ap[k as usize];
                tmp2 += ap[k as usize].conj() * x[ix as usize];
                ix += incx;
                iy += incy;
                k += 1;
            }
            y[jy as usize] += tmp * ap[(kk + j - 1) as usize].re + alpha * tmp2;
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
            let mut tmp2: Complex<T> = Complex::zero();
            y[jy as usize] += tmp * ap[kk as usize].re;
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk + 1;
            while k <= kk + n - j {
                ix += incx;
                iy += incy;
                y[iy as usize] += tmp * ap[k as usize];
                tmp2 += ap[k as usize].conj() * x[ix as usize];
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

/// HPR    performs the hermitian rank 1 operation
/// A := alpha*x*x**H + A,
/// where alpha is a real scalar, x is an n element vector and A is an n by n hermitian matrix, supplied in packed form.
/// This is [CHPR](http://www.netlib.org/lapack/explore-html/db/dcd/chpr_8f.html) and [ZHPR](http://www.netlib.org/lapack/explore-html/de/de1/zhpr_8f.html) comined in one function
#[inline]
pub fn hpr<T: Float + NumAssignOps>(
    uplo: char,
    n: isize,
    alpha: T,
    x: &[Complex<T>],
    incx: isize,
    ap: &mut [Complex<T>],
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
        panic!("hpr {}", info);
    }

    if n == 0 || alpha.is_zero() {
        return;
    }
    let mut kx = 0;
    if incx < 0 {
        kx = (-n * incx) + incx
    };
    let mut kk = 0;
    if uplo == 'u' || uplo == 'U' {
        let mut jx = kx;
        let mut j = 1;
        while j <= n {
            let x1 = x[jx as usize];
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                let mut ix = kx;
                let mut k = kk;
                while k < kk + j - 1 {
                    ap[k as usize] += tmp * x[ix as usize];
                    ix += incx;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                ap[(kk + j - 1) as usize] = Complex {
                    re: ap[(kk + j - 1) as usize].re + (x1 * tmp).re,
                    im: T::zero(),
                };
            } else {
                ap[(kk + j - 1) as usize].im = T::zero();
            }
            jx += incx;
            kk += j;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 1;
        while j <= n {
            let x1 = x[jx as usize];
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                ap[kk as usize] = Complex {
                    re: ap[kk as usize].re + (x1 * tmp).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut k = kk + 1;
                while k <= kk + n - j {
                    ix += incx;
                    ap[k as usize] += tmp * x[ix as usize];
                    k += 1;
                }
            } else {
                ap[kk as usize].im = T::zero();
            }
            jx += incx;
            kk += n - j + 1;
            j += 1;
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
    n: isize,
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
    } else if n < 0 {
        info = 2;
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
            let x1 = x[jx as usize];
            let y1 = y[jy as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                let mut ix = kx;
                let mut iy = ky;
                let mut k = kk;
                while k < kk + j - 1 {
                    ap[k as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    ix += incx;
                    iy += incy;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                ap[(kk + j - 1) as usize] = Complex {
                    re: ap[(kk + j - 1) as usize].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
            } else {
                ap[(kk + j - 1) as usize].im = T::zero();
            }
            jx += incx;
            jy += incy;
            kk += j;
            j += 1;
        }
    } else {
        let mut j = 1;
        while j <= n {
            let x1 = x[jx as usize];
            let y1 = y[jy as usize];
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                ap[kk as usize] = Complex {
                    re: ap[kk as usize].re + (x1 * tmp + y1 * tmp2).re,
                    im: T::zero(),
                };
                let mut ix = jx;
                let mut iy = jy;
                let mut k = kk + 1;
                while k <= kk + n - j {
                    ix += incx;
                    iy += incy;
                    ap[k as usize] += x[ix as usize] * tmp + y[iy as usize] * tmp2;
                    k += 1;
                }
            } else {
                ap[kk as usize].im = T::zero();
            }
            jx += incx;
            jy += incy;
            kk += n - j + 1;
            j += 1;
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
    n: isize,
    k: isize,
    a: &[Complex<T>],
    lda: isize,
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
    let noconj = trans == 't' || trans == 'T';
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
                        x[jx as usize] = tmp * a[(coords + k) as usize]
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
                let coor_aj = j * lda;
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let l = 1 - j;
                    let mut i = min(n - 1, j + k);
                    while i > j {
                        x[ix as usize] += tmp * a[((coor_aj + l + i) - 1) as usize];
                        ix -= incx;
                        i -= 1;
                    }
                    if nounit {
                        x[jx as usize] = tmp * a[coor_aj as usize];
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
            let mut tmp = x[jx as usize];
            kx -= incx;
            let mut ix = kx;
            let l = k - j;
            let coor_aj = j * lda;
            let extr_i = max(0, j - k); // evaluate once!
            if nounit {
                if noconj {
                    tmp *= a[(coor_aj + k) as usize];
                } else {
                    tmp *= a[(coor_aj + k) as usize].conj();
                }
            }
            let mut i = j - 1;
            while i >= extr_i {
                if noconj {
                    tmp += a[(coor_aj + l + i) as usize] * x[ix as usize];
                } else {
                    tmp += a[(coor_aj + l + i) as usize].conj() * x[ix as usize];
                }
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
            let mut tmp = x[jx as usize];
            kx += incx;
            let mut ix = kx;
            let l = 1 - j;
            let coor_aj = j * lda;
            let extr_i = min(n - 1, j + k); // evaluate once!
            if nounit {
                if noconj {
                    tmp *= a[coor_aj as usize];
                } else {
                    tmp *= a[coor_aj as usize].conj();
                }
            }
            let mut i = j + 1;
            while i <= extr_i {
                if noconj {
                    tmp += a[(coor_aj + l + i - 1) as usize] * x[ix as usize];
                } else {
                    tmp += a[(coor_aj + l + i - 1) as usize].conj() * x[ix as usize];
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
    n: isize,
    k: isize,
    a: &[Complex<T>],
    lda: isize,
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

    let noconj = trans == 't' || trans == 'T';
    let nounit = diag == 'n' || diag == 'N';
    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n - 1;
            while j >= 0 {
                let coords = j * lda;
                kx -= incx;
                let extri = max(0, j - k);
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let l = k - j;
                    if nounit {
                        x[jx as usize] /= a[(k + coords) as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut i = j - 1;
                    while i >= extri {
                        x[ix as usize] -= tmp * a[(coords + i + l) as usize];
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
                let extr_i = min(n - 1, j + k);
                if !x[jx as usize].is_zero() {
                    let mut ix = kx;
                    let l = 1 - j;
                    let coords = j * lda;
                    if nounit {
                        x[jx as usize] /= a[coords as usize];
                    }
                    let tmp = x[jx as usize];
                    let mut i = j + 1;
                    while i <= extr_i {
                        x[ix as usize] -= tmp * a[((coords + l + i) - 1) as usize]; // Pretty hacky maybe there is a better way?
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
            let extr_i = max(0, j - k);
            let coords = j * lda;
            let mut i = extr_i;
            while i < j {
                if noconj {
                    tmp -= a[(coords + l + i) as usize] * x[ix as usize];
                } else {
                    tmp -= a[(coords + l + i) as usize].conj() * x[ix as usize];
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[(coords + k) as usize];
                } else {
                    tmp /= a[(coords + k) as usize].conj();
                };
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
            let extr_i = min(n - 1, j + k);
            let coords = j * lda;
            let mut i = extr_i;
            while i > j {
                if noconj {
                    tmp -= a[((coords + l + i) - 1) as usize] * x[ix as usize];
                } else {
                    tmp -= a[((coords + l + i) - 1) as usize].conj() * x[ix as usize];
                };
                ix -= incx;
                i -= 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[coords as usize];
                } else {
                    tmp /= a[coords as usize].conj();
                };
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
///  x := A*x,   or   x := A**T*x,   or   x := A**H*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// This is [CTPMV](http://www.netlib.org/lapack/explore-html/d4/dbb/ctpmv_8f.html) and [ZTPMV](http://www.netlib.org/lapack/explore-html/d2/d9e/ztpmv_8f.html) comined in one function
#[inline]
pub fn tpmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
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

    let noconj = trans == 't' || trans == 'T';
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let mut k = kk;
                    while k <= kk + j - 2 {
                        x[ix as usize] += tmp * ap[k as usize];
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        x[jx as usize] = tmp * ap[(kk + j - 1) as usize];
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let mut k = kk;
                    while k >= kk - (n - (j + 1)) {
                        x[ix as usize] += tmp * ap[k as usize];
                        ix -= incx;
                        k -= 1;
                    }
                    if nounit {
                        x[jx as usize] = tmp * ap[(kk - n + j) as usize];
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
                if noconj {
                    tmp *= ap[kk as usize];
                } else {
                    tmp *= ap[kk as usize].conj();
                };
            }
            let mut k = kk - 1;
            while k > kk - j {
                ix -= incx;

                if noconj {
                    tmp += x[ix as usize] * ap[k as usize];
                } else {
                    tmp += x[ix as usize] * ap[k as usize].conj();
                };

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
                if noconj {
                    tmp *= ap[kk as usize];
                } else {
                    tmp *= ap[kk as usize].conj();
                };
            }
            let mut k = kk + 1;
            while k <= kk + n - j {
                ix += incx;
                if noconj {
                    tmp += x[ix as usize] * ap[k as usize];
                } else {
                    tmp += x[ix as usize] * ap[k as usize].conj();
                };

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
/// A*x = b,   or   A**T*x = b,   or   A**H*x = b,
/// where b and x are n element vectors and A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.
/// This is [CTPSV](http://www.netlib.org/lapack/explore-html/d8/d56/ctpsv_8f.html) and [ZTPSV](http://www.netlib.org/lapack/explore-html/da/d57/ztpsv_8f.html) comined in one function
#[inline]
pub fn tpsv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
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
    let noconj = trans == 't' || trans == 'T';
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    if nounit {
                        x[jx as usize] = tmp / ap[kk as usize];
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    if nounit {
                        x[jx as usize] = tmp / ap[kk as usize];
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
                if noconj {
                    tmp -= ap[k as usize] * x[ix as usize];
                } else {
                    tmp -= ap[k as usize].conj() * x[ix as usize];
                };
                ix += incx;
                k += 1;
            }
            if nounit {
                if noconj {
                    tmp /= ap[k as usize];
                } else {
                    tmp /= ap[k as usize].conj();
                };
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
                if noconj {
                    tmp -= ap[k as usize] * x[ix as usize];
                } else {
                    tmp -= ap[k as usize].conj() * x[ix as usize];
                };

                ix -= incx;
                k -= 1;
            }
            if nounit {
                if noconj {
                    tmp /= ap[(kk - n + j) as usize];
                } else {
                    tmp /= ap[(kk - n + j) as usize].conj();
                };
            }
            x[jx as usize] = tmp;
            jx -= incx;
            kk -= n - j + 1;
            j -= 1;
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
    n: isize,
    a: &[Complex<T>],
    lda: isize,
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

    let noconj = trans == 't' || trans == 'T';
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
                        x[jx as usize] = tmp * a[(coords + j) as usize];
                    }
                }
                j += 1;
                jx += incx;
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
                        x[jx as usize] = tmp * a[(coords + j) as usize];
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
            if noconj {
                if nounit {
                    tmp *= a[(coords + j) as usize];
                }
                while i >= 0 {
                    ix -= incx;
                    tmp += a[(coords + i) as usize] * x[ix as usize];
                    i -= 1;
                }
            } else {
                if nounit {
                    tmp *= a[(coords + j) as usize].conj();
                }
                while i >= 0 {
                    ix -= incx;
                    tmp += a[(coords + i) as usize].conj() * x[ix as usize];
                    i -= 1;
                }
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
            let mut i = j + 1;
            if noconj {
                if nounit {
                    tmp *= a[(coords + j) as usize];
                }
                while i < n {
                    ix += incx;
                    tmp += a[(coords + i) as usize] * x[ix as usize];
                    i += 1;
                }
            } else {
                if nounit {
                    tmp *= a[(coords + j) as usize].conj();
                }
                while i < n {
                    ix += incx;
                    tmp += a[(coords + i) as usize].conj() * x[ix as usize];
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
    n: isize,
    a: &[Complex<T>],
    lda: isize,
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

    let noconj = trans == 't' || trans == 'T';
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let coords = j * lda;
                    if nounit {
                        x[jx as usize] = tmp / a[(coords + j) as usize];
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
                let tmp = x[jx as usize];
                if !tmp.is_zero() {
                    let coords = j * lda;
                    if nounit {
                        x[jx as usize] = tmp / a[(coords + j) as usize];
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
            let mut ix = kx;
            let mut tmp = x[jx as usize];
            let coords = j * lda;
            let mut i = 0;
            while i < j {
                if noconj {
                    tmp -= a[(coords + i) as usize] * x[ix as usize];
                } else {
                    tmp -= a[(coords + i) as usize].conj() * x[ix as usize];
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[(coords + j) as usize];
                } else {
                    tmp /= a[(coords + j) as usize].conj();
                };
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
            let mut ix = kx;
            let mut tmp = x[jx as usize];
            let coords = j * lda;
            let mut i = n - 1;
            while i > j {
                if noconj {
                    tmp -= a[(coords + i) as usize] * x[ix as usize];
                } else {
                    tmp -= a[(coords + i) as usize].conj() * x[ix as usize];
                };
                ix -= incx;
                i -= 1;
            }
            if nounit {
                if noconj {
                    tmp /= a[(coords + j) as usize];
                } else {
                    tmp /= a[(coords + j) as usize].conj();
                };
            }
            x[jx as usize] = tmp;
            jx -= incx;
            j -= 1;
        }
    }
}
