use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};

/// AXPY constant times a vector plus a vector.
/// This is [CAXPY](http://www.netlib.org/lapack/explore-html/de/da2/caxpy_8f.html) and [ZAXPY](http://www.netlib.org/lapack/explore-html/d7/db2/zaxpy_8f.html) combined in one function.
#[inline]
pub fn axpy<T: Float + NumAssignOps>(
    n: isize,
    a: &Complex<T>,
    x: *const Complex<T>,
    incx: isize,
    y: *mut Complex<T>,
    incy: isize,
) {
    if n <= 0 {
        return;
    }
    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            *y.offset(iy) += a * *x.offset(ix);
        }
        ix += incx;
        iy += incy;
        i += 1;
    }
}

/// COPY copies a vector x to a vector y.
/// This is [CCOPY](http://www.netlib.org/lapack/explore-html/d9/dfb/ccopy_8f.html) and [ZCOPY](http://www.netlib.org/lapack/explore-html/d6/d53/zcopy_8f.html) combined in one function.
#[inline]
pub fn copy<T: Float + NumAssignOps>(
    n: isize,
    x: *const Complex<T>,
    incx: isize,
    y: *mut Complex<T>,
    incy: isize,
) {
    if n <= 0 {
        return;
    }
    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            *y.offset(iy) += *x.offset(ix);
        }
        ix += incx;
        iy += incy;
        i += 1;
    }
}

/// DOTC forms the dot product of two complex vectors DOTC = X^H * Y
/// This is [CDOTC](http://www.netlib.org/lapack/explore-html/dd/db2/cdotc_8f.html) and [ZDOTC](http://www.netlib.org/lapack/explore-html/d6/db8/zdotc_8f.html) combined in one function.
#[inline]
pub fn dotc<T: Float + NumAssignOps>(
    n: isize,
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
) -> Complex<T> {
    let mut tmp = Complex::zero();
    if n <= 0 {
        return tmp;
    }

    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            tmp += x.offset(ix).read().conj() * *y.offset(iy);
        }
        ix += incx;
        iy += incy;
        i += 1;
    }
    tmp
}

/// DOTU forms the dot product of two complex vectors DOTU = X^H * Y
/// This is [CDOTU](http://www.netlib.org/lapack/explore-html/d7/d7b/cdotu_8f.html) and [ZDOTU](http://www.netlib.org/lapack/explore-html/db/d2d/zdotu_8f.html) combined in one function.
#[inline]
pub fn dotu<T: Float + NumAssignOps>(
    n: isize,
    x: *const Complex<T>,
    incx: isize,
    y: *const Complex<T>,
    incy: isize,
) -> Complex<T> {
    let mut tmp = Complex::zero();
    if n <= 0 {
        return tmp;
    };

    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            tmp += *x.offset(ix) * *y.offset(iy);
        }
        ix += incx;
        iy += incy;
        i += 1;
    }
    tmp
}

/// ROT applies a plane rotation, where the cos and sin (c and s) are real and the vectors cx and cy are complex.
/// This is [CSROT](http://www.netlib.org/lapack/explore-html/d1/dbb/csrot_8f.html) and [ZDROT](http://www.netlib.org/lapack/explore-html/d4/de9/zdrot_8f.html) combined in one function.
#[inline]
pub fn rot<T: Float + NumAssignOps>(
    n: isize,
    x: *mut Complex<T>,
    incx: isize,
    y: *mut Complex<T>,
    incy: isize,
    c: T,
    s: T,
) {
    if n <= 0 {
        return;
    }
    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;
    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            let tmp = *x.offset(ix) * c + *y.offset(iy) * s;
            *y.offset(iy) = *y.offset(iy) * c - *x.offset(ix) * s;
            *x.offset(ix) = tmp;
        }
        i += 1;
        ix += incx;
        iy += incy;
    }
}

/// ROTG determines a complex Givens rotation.
/// This is [CROTG](http://www.netlib.org/lapack/explore-html/dc/de6/crotg_8f.html) and [ZROTG](http://www.netlib.org/lapack/explore-html/dc/dfe/zrotg_8f.html) combined in one function.
#[inline]
pub fn rotg<T: Float + NumAssignOps>(
    a: &mut Complex<T>,
    b: &mut Complex<T>,
    c: &mut T,
    s: &mut Complex<T>,
) {
    if a.is_zero() {
        *c = T::zero();
        *s = Complex::one();
        *a = *b;
        return;
    }
    let scale = a.norm() + b.norm();
    let norm = scale * ((*a / scale).norm().powi(2) + (*b / scale).norm().powi(2)).sqrt();
    let alpha = *a / a.norm();
    *c = a.norm() / norm;
    *s = alpha * b.conj() / norm;
    *a = alpha * norm;
}

/// SCAL scales a complex vector by a constant.
/// This is [CSCAL](http://www.netlib.org/lapack/explore-html/dc/d81/cscal_8f.html) and [ZSCAL](http://www.netlib.org/lapack/explore-html/d2/d74/zscal_8f.html) combined in one function.
#[inline]
pub fn scal<T: Float + NumAssignOps>(n: isize, a: Complex<T>, x: *mut Complex<T>, incx: isize) {
    if n <= 0 {
        return;
    }
    if incx <= 0 {
        return;
    }

    let mut i = 0;
    let nincx = n * incx;
    while i < nincx {
        unsafe {
            *x.offset(i) *= a;
        }
        i += incx;
    }
}

/// SSCAL scales a complex vector by a real constant.
/// This is [CSSCAL](http://www.netlib.org/lapack/explore-html/de/d5e/csscal_8f.html) and [ZDSCAL](http://www.netlib.org/lapack/explore-html/dd/d76/zdscal_8f.html) combined in one function.
#[inline]
pub fn sscal<T: Float + NumAssignOps>(n: isize, a: T, x: *mut Complex<T>, incx: isize) {
    if n <= 0 {
        return;
    }
    if incx <= 0 {
        return;
    }

    let mut i = 0;
    let nincx = n * incx;
    while i < nincx {
        unsafe {
            if !a.is_zero() {
                *x.offset(i) *= a;
            } else {
                *x.offset(i) = Complex::zero();
            }
        }
        i += incx;
    }
}

/// SWAP interchanges two vectors.
/// This is [CSWAP](http://www.netlib.org/lapack/explore-html/d1/d44/cswap_8f.html) and [ZSWAP](http://www.netlib.org/lapack/explore-html/d3/dc0/zswap_8f.html) combined in one function.
#[inline]
pub fn swap<T: Float + NumAssignOps>(
    n: isize,
    x: *mut Complex<T>,
    incx: isize,
    y: *mut Complex<T>,
    incy: isize,
) {
    if n <= 0 {
        return;
    }

    let mut ix = 0;
    let mut iy = 0;
    let mut i = 0;

    if incx < 0 {
        ix = (-n * incx) + incx;
    }
    if incy < 0 {
        iy = (-n * incy) + incy;
    }
    while i < n {
        unsafe {
            let tmp = *x.offset(ix);
            *x.offset(ix) = *y.offset(iy);
            *y.offset(iy) = tmp;
        }
        ix += incx;
        iy += incy;
        i += 1;
    }
}

/// IAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
/// This is [ICAMAX](http://www.netlib.org/lapack/explore-html/dd/d51/icamax_8f.html) and [IZAMAX](http://www.netlib.org/lapack/explore-html/d0/da5/izamax_8f.html) combined in one function.
#[inline]
pub fn iamax<T: Float + NumAssignOps>(n: isize, x: *const Complex<T>, incx: isize) -> isize {
    let mut iamax = 0;
    if n <= 0 || incx <= 0 {
        return iamax;
    }
    iamax = 1;
    if n == 1 {
        return iamax;
    }

    // FIXME use l1_norm()
    let mut max;
    unsafe { max = x.read().norm() };
    let mut i = 2;
    if incx == 1 {
        while i < n {
            // FIXME use l1_norm()
            let tmp;
            unsafe {
                tmp = x.offset(i).read().norm();
            }
            i += 1;
            if tmp > max {
                iamax = i;
                max = tmp;
            }
        }
    } else {
        let mut ix = 1 + incx;
        while i < n {
            // FIXME use l1_norm()
            let tmp;
            unsafe {
                tmp = x.offset(ix).read().norm();
            }
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

/// ASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and returns a single or double precision result.
/// This is [SCASUM](http://www.netlib.org/lapack/explore-html/db/d53/scasum_8f.html) and [DZASUM](http://www.netlib.org/lapack/explore-html/df/d0f/dzasum_8f.html) combined in one function.
#[inline]
pub fn asum<T: Float + NumAssignOps>(n: isize, x: *const Complex<T>, incx: isize) -> T {
    let mut sum = T::zero();
    if n <= 0 || incx <= 0 {
        return sum;
    }

    let mut i = 0;
    let nincx = n * incx;
    while i < nincx {
        unsafe {
            let Complex { re, im } = *x.offset(i);
            // FIXME use l1_norm()
            sum += re.abs() + im.abs();
        }
        i += incx;
    }
    sum
}

/// NRM2 returns the euclidean norm of a vector via the function name, so that NRM2 := sqrt( x**H*x )
/// This is [SCNRM2](http://www.netlib.org/lapack/explore-html/db/d66/scnrm2_8f.html) and [DZNRM2](http://www.netlib.org/lapack/explore-html/d9/d19/dznrm2_8f.html) combined in one function.
pub fn nrm2<T: Float + NumAssignOps>(n: isize, x: *const Complex<T>, incx: isize) -> T {
    if n <= 0 || incx <= 0 {
        return T::zero();
    }
    //  The following loop is equivalent to this call to the LAPACK auxiliary routine:
    //  CALL CLASSQ( N, X, INCX, SCALE, SSQ )
    let mut scale = T::zero();
    let mut ssq = T::one();
    let mut ix = 0;
    while ix < n * incx {
        let Complex { re, im } = unsafe { *x.offset(ix) };
        if !re.is_zero() {
            let tmp = re.abs();
            if scale < tmp {
                let t1 = scale / tmp;
                ssq = T::one() + ssq * t1 * t1;
                scale = tmp;
            } else {
                let t1 = tmp / scale;
                ssq += t1 * t1;
            }
        }
        if !im.is_zero() {
            let tmp = im.abs();
            if scale < tmp {
                let t1 = scale / tmp;
                ssq = T::one() + ssq * t1 * t1;
                scale = tmp;
            } else {
                let t1 = tmp / scale;
                ssq += t1 * t1;
            }
        }
        ix += incx
    }
    scale * ssq.sqrt()
}
