use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::{max, min};

fn multiply<T: Float + NumAssignOps>(
    left: *mut Complex<T>,
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
            unsafe {
                *left.offset(index) *= right;
            }
            index += inc;
            i += 1;
        }
    }
}

fn zero<T: Float + NumAssignOps>(left: *mut Complex<T>, len: isize, mut index: isize, inc: isize) {
    let mut i = 0;
    while i < len {
        unsafe {
            *left.offset(index) = Complex::zero();
        }
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
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    x: *const Complex<T>,
    incx: isize,
    beta: Complex<T>,
    y: *mut Complex<T>,
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

    let noconj = trans == 't' || trans == 'T';

    let mut kx = 0;
    let mut ky = 0;
    if incx < 0 {
        kx = (-lenx * incx) + incx
    };
    if incy < 0 {
        ky = (-leny * incy) + incy
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
            let tmp = unsafe { *x.offset(jx) * alpha };
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
            let k = ku - j; // apparently we want this to be negative sometimes.
            let aj = j * lda;
            let mut i = max(0, j - ku);
            while i < min(m, j + kl + 1) {
                if noconj {
                    unsafe {
                        tmp += *a.offset(aj + k + i) * *x.offset(ix);
                    }
                } else {
                    unsafe {
                        tmp += a.offset(aj + k + i).read().conj() * *x.offset(ix);
                    }
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
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    x: *const Complex<T>,
    incx: isize,
    beta: Complex<T>,
    y: *mut Complex<T>,
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
        kx = (-lenx * incx) + incx
    };
    if incy < 0 {
        ky = (-leny * incx) + incy
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
            let mut tmp: Complex<T> = Complex::zero();
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                if noconj {
                    unsafe {
                        tmp += *a.offset(aj + i) * *x.offset(ix);
                    }
                } else {
                    unsafe {
                        tmp += a.offset(aj + i).read().conj() * *x.offset(ix);
                    }
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

/// GERC  performs the rank 1 operation
/// A := alpha*x*y**H + A,
/// where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n matrix.
/// This is [CGERC](http://www.netlib.org/lapack/explore-html/dd/d84/cgerc_8f.html) and [ZGERC](http://www.netlib.org/lapack/explore-html/d3/dad/zgerc_8f.html) comined in one function
#[inline]
pub fn gerc<T: Float + NumAssignOps>(
    m: isize,
    n: isize,
    alpha: Complex<T>,
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
    a: *mut Complex<T>,
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
        panic!("gerc {}", info);
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
        let tmp = unsafe { *y.offset(jy) };
        if !tmp.is_zero() {
            let tmp = alpha * tmp.conj();
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp };
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
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
    a: *mut Complex<T>,
    lda: isize,
) {
    let mut info = 0;
    if m < 0 {
        info = 5;
    } else if n < 0 {
        info = 5;
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
        let tmp = unsafe { *y.offset(jy) };
        if !tmp.is_zero() {
            let tmp = alpha * tmp;
            let mut ix = kx;
            let aj = j * lda;
            let mut i = 0;
            while i < m {
                unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp };
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
    a: *const Complex<T>,
    lda: isize,
    x: *const Complex<T>,
    incx: isize,
    beta: Complex<T>,
    y: *mut Complex<T>,
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
        panic!("hbmv {}", info);
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
    if uplo == 'u' || uplo == 'U' {
        let mut j = 0;
        while j < n {
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = max(0, j - k);
            while i < j {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + k - j + i);
                    tmp2 += a.offset(aj + k - j + i).read().conj() * *x.offset(ix);
                }
                ix += incx;
                iy += incy;
                i += 1;
            }
            unsafe { *y.offset(iy) += tmp * a.offset(aj + k).read().re + (alpha * tmp2) };
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
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2: Complex<T> = Complex::zero();
            let aj = j * lda;
            unsafe { *y.offset(jy) += tmp * a.offset(aj).read().re };
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj - j + i);
                    tmp2 += a.offset(aj - j + i).read().conj() * *x.offset(ix);
                }
                i += 1;
            }
            unsafe { *y.offset(jy) += alpha * tmp2 };
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
    a: *const Complex<T>,
    lda: isize,
    x: *const Complex<T>,
    incx: isize,
    beta: Complex<T>,
    y: *mut Complex<T>,
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
        panic!("hemv {}", info);
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
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + i);
                    tmp2 += a.offset(aj + i).read().conj() * *x.offset(ix);
                }
                ix += incx;
                iy += incy;
                i += 1;
            }
            unsafe { *y.offset(jy) += tmp * a.offset(aj + j).read().re + (alpha * tmp2) };
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        while j < n {
            let aj = j * lda;
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2: Complex<T> = Complex::zero();
            unsafe { *y.offset(jy) += tmp * a.offset(aj + j).read().re };
            let mut ix = jx;
            let mut iy = jy;
            let mut i = j + 1;
            while i < n {
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *a.offset(aj + i);
                    tmp2 += a.offset(aj + i).read().conj() * *x.offset(ix);
                }
                i += 1;
            }
            unsafe { *y.offset(jy) += alpha * tmp2 };
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
    x: *const Complex<T>,
    incx: isize,
    a: *mut Complex<T>,
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
            let aj = j * lda;
            let tmp = unsafe { *x.offset(jx) };
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let mut ix = kx;
                let mut i = 0;
                while i < j {
                    unsafe {
                        *a.offset(aj + i) += tmp * *x.offset(ix);
                    }
                    ix += incx;
                    i += 1;
                }
                let tmp = unsafe { tmp * *x.offset(jx) };
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re + tmp.re,
                        im: T::zero(),
                    };
                }
            } else {
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re,
                        im: T::zero(),
                    };
                }
            }
            jx += incx;
            j += 1;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let tmp = unsafe { *x.offset(jx) };
            if !tmp.is_zero() {
                let tmp = tmp.conj() * alpha;
                let tmp2 = unsafe { tmp * *x.offset(jx) };
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re + tmp2.re,
                        im: T::zero(),
                    };
                }
                let mut ix = jx;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    unsafe {
                        *a.offset(aj + i) += tmp * *x.offset(ix);
                    }
                    i += 1;
                }
            } else {
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re,
                        im: T::zero(),
                    };
                }
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
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
    a: *mut Complex<T>,
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
            let aj = j * lda;
            let y1 = unsafe { *y.offset(jy) };
            let x1 = unsafe { *x.offset(jx) };
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();

                let mut ix = kx;
                let mut iy = ky;
                let mut i = 0;
                while i < j {
                    unsafe {
                        *a.offset(aj + i) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2;
                    }
                    ix += incx;
                    iy += incy;
                    i += 1;
                }
                // FIXME
                // maybe only calculate the real.
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re + (x1 * tmp + y1 * tmp2).re,
                        im: T::zero(),
                    };
                }
            } else {
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re,
                        im: T::zero(),
                    };
                }
            }
            jx += incx;
            jy += incy;
            j += 1;
        }
    } else {
        while j < n {
            let aj = j * lda;
            let y1 = unsafe { *y.offset(jy) };
            let x1 = unsafe { *x.offset(jx) };
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re + (x1 * tmp + y1 * tmp2).re,
                        im: T::zero(),
                    };
                }
                let mut ix = jx;
                let mut iy = jy;
                let mut i = j + 1;
                while i < n {
                    ix += incx;
                    iy += incy;
                    unsafe { *a.offset(aj + i) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                    i += 1;
                }
            } else {
                unsafe {
                    *a.offset(aj + j) = Complex {
                        re: a.offset(aj + j).read().re,
                        im: T::zero(),
                    };
                }
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
    ap: *const Complex<T>,
    x: *const Complex<T>,
    incx: isize,
    beta: Complex<T>,
    y: *mut Complex<T>,
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
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = kx;
            let mut iy = ky;
            let mut k = kk;
            kk += j;
            while k < kk - 1 {
                // KK + J - 2
                unsafe {
                    *y.offset(iy) += tmp * *ap.offset(k);
                    tmp2 += ap.offset(k).read().conj() * *x.offset(ix);
                };
                ix += incx;
                iy += incy;
                k += 1;
            }
            unsafe { *y.offset(jy) += tmp * ap.offset(kk - 1).read().re + alpha * tmp2 };
            jx += incx;
            jy += incy;
        }
    } else {
        let mut jx = kx;
        let mut jy = ky;
        let mut j = 0;
        while j < n {
            j += 1;
            let tmp = unsafe { alpha * *x.offset(jx) };
            let mut tmp2: Complex<T> = Complex::zero();
            let mut ix = jx;
            let mut iy = jy;
            let mut k = kk;
            kk += n - j;
            unsafe { *y.offset(jy) += tmp * ap.offset(k).read().re };
            while k < kk {
                k += 1;
                ix += incx;
                iy += incy;
                unsafe {
                    *y.offset(iy) += tmp * *ap.offset(k);
                    tmp2 += ap.offset(k).read().conj() * *x.offset(ix);
                };
            }
            unsafe { *y.offset(jy) += alpha * tmp2 };
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
    n: isize,
    alpha: T,
    x: *const Complex<T>,
    incx: isize,
    ap: *mut Complex<T>,
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
        let mut j = 0;
        while j < n {
            let x1 = unsafe { *x.offset(jx) };
            let mut k = kk;
            kk += j;
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                let mut ix = kx;
                while k < kk {
                    unsafe { *ap.offset(k) += tmp * *x.offset(ix) };
                    ix += incx;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                unsafe {
                    *ap.offset(kk) = Complex {
                        re: ap.offset(kk).read().re + (x1 * tmp).re,
                        im: T::zero(),
                    };
                }
            } else {
                unsafe {
                    *ap.offset(kk) = Complex {
                        re: ap.offset(kk).read().re,
                        im: T::zero(),
                    };
                };
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
            let x1 = unsafe { *x.offset(jx) };
            let mut k = kk;
            kk += n - j;
            if !x1.is_zero() {
                let tmp = x1.conj() * alpha;
                unsafe {
                    *ap.offset(k) = Complex {
                        re: ap.offset(k).read().re + (x1 * tmp).re,
                        im: T::zero(),
                    };
                }
                let mut ix = jx;
                while k < kk {
                    k += 1;
                    ix += incx;
                    unsafe { *ap.offset(k) += tmp * *x.offset(ix) };
                }
            } else {
                unsafe {
                    *ap.offset(k) = Complex {
                        re: ap.offset(k).read().re,
                        im: T::zero(),
                    };
                };
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
    n: isize,
    alpha: Complex<T>,
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
    ap: *mut Complex<T>,
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
        let mut j = 0;
        while j < n {
            let x1 = unsafe { *x.offset(jx) };
            let y1 = unsafe { *y.offset(jy) };
            let mut k = kk;
            kk += j;
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                let mut ix = kx;
                let mut iy = ky;
                while k < kk {
                    unsafe { *ap.offset(k) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                    ix += incx;
                    iy += incy;
                    k += 1;
                }
                // FIXME
                // maybe only calculate the real.
                unsafe {
                    *ap.offset(kk) = Complex {
                        re: ap.offset(kk).read().re + (x1 * tmp + y1 * tmp2).re,
                        im: T::zero(),
                    };
                }
            } else {
                unsafe {
                    *ap.offset(kk) = Complex {
                        re: ap.offset(kk).read().re,
                        im: T::zero(),
                    };
                };
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
            let x1 = unsafe { *x.offset(jx) };
            let y1 = unsafe { *y.offset(jy) };
            let mut k = kk;
            kk += n - j;
            if !x1.is_zero() || !y1.is_zero() {
                let tmp = alpha * y1.conj();
                let tmp2 = (alpha * x1).conj();
                // FIXME
                // maybe only calculate the real.
                unsafe {
                    *ap.offset(k) = Complex {
                        re: ap.offset(k).read().re + (x1 * tmp + y1 * tmp2).re,
                        im: T::zero(),
                    };
                }
                let mut ix = jx;
                let mut iy = jy;
                // FIXME Maybe we could move kk += n - j + 1 up here
                while k < kk {
                    k += 1;
                    ix += incx;
                    iy += incy;
                    unsafe { *ap.offset(k) += *x.offset(ix) * tmp + *y.offset(iy) * tmp2 };
                }
            } else {
                unsafe {
                    *ap.offset(k) = Complex {
                        re: ap.offset(k).read().re,
                        im: T::zero(),
                    };
                };
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
    n: isize,
    k: isize,
    a: *const Complex<T>,
    lda: isize,
    x: *mut Complex<T>,
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
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let aj = j * lda;
                    let mut i = max(0, j - k);
                    while i < j {
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj + i + k - j) };
                        ix += incx;
                        i += 1;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) = tmp * *a.offset(aj + k) }
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
            kx += (n - 1) * incx;
            let mut jx = kx;
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let aj = j * lda;
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let mut ix = kx;
                    let mut i = min(n - 1, j + k); // NOTE
                    while i >= j + 1 {
                        unsafe { *x.offset(ix) += tmp * *a.offset(aj - j + i) };
                        ix -= incx;
                        i -= 1;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) = tmp * *a.offset(aj) };
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
            let mut tmp = unsafe { *x.offset(jx) };
            kx -= incx;
            let mut ix = kx;
            let aj = j * lda;
            if nounit {
                if noconj {
                    unsafe { tmp *= *a.offset(aj + k) };
                } else {
                    unsafe { tmp *= a.offset(aj + k).read().conj() };
                }
            }
            let mut i = j;
            // FIXME figure out if there is a better way
            while i >= max(1, j - k + 1) {
                i -= 1;
                if noconj {
                    unsafe { tmp += *a.offset(aj + i + k - j) * *x.offset(ix) };
                } else {
                    unsafe { tmp += a.offset(aj + i + k - j).read().conj() * *x.offset(ix) };
                }
                ix -= incx;
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
        }
    } else {
        let mut jx = kx;
        let mut j = 0;
        while j < n {
            let mut tmp = unsafe { *x.offset(jx) };
            kx += incx;
            let aj = j * lda;
            if nounit {
                if noconj {
                    unsafe { tmp *= *a.offset(aj) };
                } else {
                    unsafe { tmp *= a.offset(aj).read().conj() };
                }
            }
            let mut ix = kx;
            let mut i = j + 1;
            while i < min(n, j + k + 1) {
                if noconj {
                    unsafe { tmp += *a.offset(aj - j + i) * *x.offset(ix) };
                } else {
                    unsafe { tmp += a.offset(aj - j + i).read().conj() * *x.offset(ix) };
                }

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
    a: *const Complex<T>,
    lda: isize,
    x: *mut Complex<T>,
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
            let mut j = n;
            while j >= 1 {
                j -= 1;
                kx -= incx;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let mut ix = kx;
                    let aj = j * lda;
                    if nounit {
                        unsafe { *x.offset(jx) /= *a.offset(aj + k) };
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut i = j;
                    while i > max(0, j - k) {
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
                    if nounit {
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
                if noconj {
                    unsafe { tmp -= *a.offset(aj + k - j + i) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= a.offset(aj + k - j + i).read().conj() * *x.offset(ix) };
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    unsafe { tmp /= *a.offset(aj + k) };
                } else {
                    unsafe { tmp /= a.offset(aj + k).read().conj() };
                };
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
            let aj = j * lda;
            let mut tmp = unsafe { *x.offset(jx) };
            let mut ix = kx;
            //FIXME Maybe we could do this in a diffrent way
            let mut i = min(n - 1, j + k);
            while i > j {
                if noconj {
                    unsafe { tmp -= *a.offset(aj - j + i) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= a.offset(aj - j + i).read().conj() * *x.offset(ix) };
                };
                i -= 1;
                ix -= incx;
            }
            //
            if nounit {
                if noconj {
                    unsafe { tmp /= *a.offset(aj) };
                } else {
                    unsafe { tmp /= a.offset(aj).read().conj() };
                };
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
///  x := A*x,   or   x := A**T*x,   or   x := A**H*x,
/// where x is an n element vector and  A is an n by n unit, or non-unit, upper or lower triangular matrix, supplied in packed form.
/// This is [CTPMV](http://www.netlib.org/lapack/explore-html/d4/dbb/ctpmv_8f.html) and [ZTPMV](http://www.netlib.org/lapack/explore-html/d2/d9e/ztpmv_8f.html) comined in one function
#[inline]
pub fn tpmv<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    diag: char,
    n: isize,
    ap: *const Complex<T>,
    x: *mut Complex<T>,
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
            let mut j = 0;
            while j < n {
                j += 1;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = kx;
                    let mut k = kk;
                    while k < kk + j - 1 {
                        unsafe { *x.offset(ix) += tmp * *ap.offset(k) };
                        ix += incx;
                        k += 1;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) *= *ap.offset(kk + j - 1) };
                    }
                }
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
                //FIXME figure out a better way
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = kx;
                    let mut k = kk - 1;
                    while k > kk - (n - j) {
                        unsafe { *x.offset(ix) += tmp * *ap.offset(k) };
                        ix -= incx;
                        k -= 1;
                    }
                    if nounit {
                        unsafe { *x.offset(jx) *= *ap.offset(kk - n + j) };
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
                if noconj {
                    tmp *= unsafe { *ap.offset(k) };
                } else {
                    tmp *= unsafe { ap.offset(k).read().conj() };
                };
            }
            while k > kk {
                ix -= incx;
                k -= 1;
                if noconj {
                    unsafe { tmp += *ap.offset(k) * *x.offset(ix) };
                } else {
                    unsafe { tmp += ap.offset(k).read().conj() * *x.offset(ix) };
                };
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
                if noconj {
                    unsafe { tmp *= *ap.offset(k) };
                } else {
                    unsafe { tmp *= ap.offset(k).read().conj() };
                };
            }
            while k < kk {
                ix += incx;
                k += 1;
                if noconj {
                    unsafe { tmp += *ap.offset(k) * *x.offset(ix) };
                } else {
                    unsafe { tmp += ap.offset(k).read().conj() * *x.offset(ix) };
                };
            }
            unsafe { *x.offset(jx) = tmp };
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
    n: isize,
    ap: *const Complex<T>,
    x: *mut Complex<T>,
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
            let mut kk = n * (n + 1) / 2;
            let mut jx = kx + (n - 1) * incx;
            let mut j = n;
            while j >= 1 {
                kk -= j;
                j -= 1;
                if unsafe { !x.offset(jx).read().is_zero() } {
                    let mut k = kk + j;
                    if nounit {
                        unsafe { *x.offset(jx) /= *ap.offset(k) };
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
                        unsafe { *x.offset(jx) /= *ap.offset(k) };
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
            kk += j;
            while k < kk - 1 {
                if noconj {
                    unsafe { tmp -= *ap.offset(k) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= ap.offset(k).read().conj() * *x.offset(ix) };
                };
                ix += incx;
                k += 1;
            }
            if nounit {
                // NOTE k == kk + j - 1
                if noconj {
                    unsafe { tmp /= *ap.offset(k) };
                } else {
                    unsafe { tmp /= ap.offset(k).read().conj() };
                };
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
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
                if noconj {
                    unsafe { tmp -= *ap.offset(k) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= ap.offset(k).read().conj() * *x.offset(ix) };
                };

                ix -= incx;
                k -= 1;
            }
            if nounit {
                // NOTE k == kk - n + j
                if noconj {
                    unsafe { tmp /= *ap.offset(k) };
                } else {
                    unsafe { tmp /= ap.offset(k).read().conj() };
                };
            }
            unsafe { *x.offset(jx) = tmp };
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
    n: isize,
    a: *const Complex<T>,
    lda: isize,
    x: *mut Complex<T>,
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
                        unsafe { *x.offset(jx) = tmp * *a.offset(aj + j) };
                    }
                }
                j += 1;
                jx += incx;
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
                        unsafe { *x.offset(jx) = tmp * *a.offset(aj + j) };
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
            if noconj {
                if nounit {
                    unsafe { tmp *= *a.offset(aj + j) };
                }
                while i >= 1 {
                    i -= 1;
                    ix -= incx;
                    unsafe { tmp += *a.offset(aj + i) * *x.offset(ix) };
                }
            } else {
                if nounit {
                    unsafe { tmp *= a.offset(aj + j).read().conj() };
                }
                while i >= 1 {
                    i -= 1;
                    ix -= incx;
                    unsafe { tmp += a.offset(aj + i).read().conj() * *x.offset(ix) };
                }
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
            let mut i = j + 1;
            if noconj {
                if nounit {
                    unsafe { tmp *= *a.offset(aj + j) };
                }
                while i < n {
                    ix += incx;
                    unsafe { tmp += *a.offset(aj + i) * *x.offset(ix) };
                    i += 1;
                }
            } else {
                if nounit {
                    unsafe { tmp *= a.offset(aj + j).read().conj() };
                }
                while i < n {
                    ix += incx;
                    unsafe { tmp += a.offset(aj + i).read().conj() * *x.offset(ix) };
                    i += 1;
                }
            }
            unsafe { *x.offset(jx) = tmp };
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
    a: *const Complex<T>,
    lda: isize,
    x: *mut Complex<T>,
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
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let aj = j * lda;
                    if nounit {
                        unsafe { *x.offset(jx) = tmp / *a.offset(aj + j) };
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
                let tmp = unsafe { *x.offset(jx) };
                if !tmp.is_zero() {
                    let aj = j * lda;
                    if nounit {
                        unsafe { *x.offset(jx) = tmp / *a.offset(aj + j) };
                    }
                    let tmp = unsafe { *x.offset(jx) };
                    let mut ix = jx;
                    let mut i = j + 1;
                    while i < n {
                        ix += incx;
                        unsafe { *x.offset(ix) -= tmp * *a.offset(aj + i) };
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
            let mut tmp = unsafe { *x.offset(jx) };
            let aj = j * lda;
            let mut i = 0;
            while i < j {
                if noconj {
                    unsafe { tmp -= *a.offset(aj + i) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= a.offset(aj + i).read().conj() * *x.offset(ix) };
                };
                ix += incx;
                i += 1;
            }
            if nounit {
                if noconj {
                    unsafe { tmp /= *a.offset(aj + j) };
                } else {
                    unsafe { tmp /= a.offset(aj + j).read().conj() };
                };
            }
            unsafe { *x.offset(jx) = tmp };
            jx += incx;
            j += 1;
        }
    } else {
        kx += (n - 1) * incx;
        let mut jx = kx;
        let mut j = n;
        while j >= 1 {
            j -= 1;
            let mut ix = kx;
            let mut tmp = unsafe { *x.offset(jx) };
            let aj = j * lda;
            let mut i = n;
            while i > j + 1 {
                i -= 1;
                if noconj {
                    unsafe { tmp -= *a.offset(aj + i) * *x.offset(ix) };
                } else {
                    unsafe { tmp -= a.offset(aj + i).read().conj() * *x.offset(ix) };
                };
                ix -= incx;
            }
            if nounit {
                if noconj {
                    unsafe { tmp /= *a.offset(aj + j) };
                } else {
                    unsafe { tmp /= a.offset(aj + j).read().conj() };
                };
            }
            unsafe { *x.offset(jx) = tmp };
            jx -= incx;
        }
    }
}
