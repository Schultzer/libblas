use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::max;

fn multiply<T: Float + NumAssignOps>(
    left: *mut Complex<T>,
    right: Complex<T>,
    ld: isize,
    n: isize,
    m: isize,
) {
    if right.is_zero() {
        zero(left, ld, n, m);
    } else {
        let mut j = 0;
        while j < n {
            let coords = j * ld;
            let mut i = 0;
            while i < m {
                unsafe { *left.offset(coords + i) *= right };
                i += 1;
            }
            j += 1;
        }
    }
}

fn zero<T: Float + NumAssignOps>(left: *mut Complex<T>, ld: isize, n: isize, m: isize) {
    let mut j = 0;
    while j < n {
        let coords = j * ld;
        let mut i = 0;
        while i < m {
            unsafe { *left.offset(coords + i) = Complex::zero() };
            i += 1;
        }
        j += 1;
    }
}

/// CGEMM  performs one of the matrix-matrix operations
/// C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
/// alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
/// This is [CGEMM](http://www.netlib.org/lapack/explore-html/d6/d5b/cgemm_8f.html) and [ZGEMM](http://www.netlib.org/lapack/explore-html/d7/d76/zgemm_8f.html) comined in one function
#[inline]
pub fn gemm<T: Float + NumAssignOps>(
    trans_a: char,
    trans_b: char,
    m: isize,
    n: isize,
    k: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *const Complex<T>,
    ldb: isize,
    beta: Complex<T>,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let nrowa = if trans_a == 'n' || trans_a == 'N' {
        m
    } else {
        k
    };
    let nrowb = if trans_b == 'n' || trans_b == 'N' {
        k
    } else {
        n
    };

    let mut info = 0;
    if !(trans_a == 'T' || trans_a == 't' || trans_a == 'n' || trans_a == 'C' || trans_a == 'c') {
        info = 1;
    } else if !(trans_b == 'T'
        || trans_b == 't'
        || trans_b == 'n'
        || trans_b == 'C'
        || trans_b == 'c')
    {
        info = 2;
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if k < 0 {
        info = 5;
    } else if lda < max(1, nrowa) {
        info = 8;
    } else if ldb < max(1, nrowb) {
        info = 10;
    } else if ldc < max(1, m) {
        info = 13;
    }
    if info != 0 {
        panic!("gemm {}", info);
    }

    let alpha_is_zero = alpha.is_zero();
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();

    if m == 0 || n == 0 || (alpha_is_zero || k == 0) && beta_is_one {
        return;
    }

    if alpha_is_zero {
        multiply(c, beta, ldc, n, m);
        return;
    }

    if (trans_a == 'n' || trans_a == 'N') && (trans_b == 'n' || trans_b == 'N') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let bj = j * ldb;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) *= beta };
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let tmp = unsafe { alpha * *b.offset(bj + l) };
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *a.offset(al + i) };
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'n' || trans_a == 'N') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) *= beta };
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let bl = l * ldb;
                let al = l * lda;
                let tmp = unsafe { alpha * b.offset(bl + j).read().conj() };
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *a.offset(al + i) };
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'n' || trans_a == 'N') && (trans_b == 't' || trans_b == 'T') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) *= beta };
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let bl = l * ldb;
                let al = l * lda;
                let tmp = unsafe { alpha * *b.offset(bl + j) };
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *a.offset(al + i) };
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'c' || trans_a == 'C') && (trans_b == 'n' || trans_b == 'N') {
        let mut j = 0;
        while j < n {
            let bj = j * ldb;
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe { tmp += a.offset(ai + l).read().conj() * *b.offset(bj + l) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) }
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'c' || trans_a == 'C') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let bl = l * ldb;
                    unsafe {
                        tmp += a.offset(ai + l).read().conj() * b.offset(bl + j).read().conj()
                    };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'c' || trans_a == 'C') && (trans_b == 't' || trans_b == 'T') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let bl = l * ldb;
                    unsafe { tmp += a.offset(ai + l).read().conj() * *b.offset(bl + j) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 't' || trans_a == 'T') && (trans_b == 'n' || trans_b == 'N') {
        let mut j = 0;
        while j < n {
            let bj = j * ldb;
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe { tmp += *a.offset(ai + l) * *b.offset(bj + l) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 't' || trans_a == 'T') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let bl = l * ldb;
                    unsafe { tmp += *a.offset(ai + l) * b.offset(bl + j).read().conj() };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
        return;
    } else {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let mut i = 0;
            while i < m {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let bl = l * ldb;
                    unsafe { tmp += *a.offset(ai + l) * *b.offset(bl + j) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
    }
    return;
}

/// CSYMM  performs one of the matrix-matrix operations
/// C := alpha*A*B + beta*C or C := alpha*B*A + beta*C,
/// where  alpha and beta are scalars, A is a symmetric matrix and  B and C are m by n matrices.
/// This is [CSYMM](http://www.netlib.org/lapack/explore-html/db/d59/csymm_8f.html) and [ZSYMM](http://www.netlib.org/lapack/explore-html/df/d51/zsymm_8f.html) comined in one function
#[inline]
pub fn symm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *const Complex<T>,
    ldb: isize,
    beta: Complex<T>,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let nrowa = if side == 'l' || side == 'L' { m } else { n };
    let upper = uplo == 'u' || uplo == 'U';
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldb < max(1, m) {
        info = 9;
    } else if ldc < max(1, m) {
        info = 12;
    }
    if info != 0 {
        panic!("symm {}", info);
    }

    let alpha_is_zero = alpha.is_zero();
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();

    if m == 0 || n == 0 || (alpha_is_zero && beta_is_one) {
        return;
    }

    if alpha_is_zero {
        multiply(c, beta, ldc, n, m);
        return;
    }

    if side == 'l' || side == 'L' {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let cj = j * ldc;
                let mut i = 0;
                while i < m {
                    let ai = i * ldc;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        unsafe {
                            *c.offset(cj + k) += tmp * *a.offset(ai + k);
                            tmp2 += *b.offset(bj + k) * *a.offset(ai + k);
                        }
                        k += 1
                    }
                    tmp = unsafe { tmp * *a.offset(ai + i) + alpha * tmp2 };
                    if !beta_is_zero {
                        unsafe { tmp += beta * *c.offset(cj + i) };
                    }
                    unsafe { *c.offset(cj + i) = tmp };
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let cj = j * ldc;
                let mut i = m;
                while i >= 1 {
                    i -= 1;
                    let ai = i * lda;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        unsafe {
                            *c.offset(cj + k) += tmp * *a.offset(ai + k);
                            tmp2 += *b.offset(bj + k) * *a.offset(ai + k);
                        }
                        k += 1
                    }
                    tmp = unsafe { tmp * *a.offset(ai + i) + alpha * tmp2 };
                    if !beta_is_zero {
                        unsafe { tmp += beta * *c.offset(cj + i) };
                    }
                    unsafe { *c.offset(cj + i) = tmp };
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let bj = j * ldb;
            let cj = j * ldc;
            let mut tmp = unsafe { alpha * *a.offset(aj + j) };
            let mut i = 0;
            while i < m {
                let mut tmp2 = unsafe { tmp * *b.offset(bj + i) };
                if !beta_is_zero {
                    unsafe { tmp2 += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp2 };
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let ak = k * lda;
                let bk = k * ldb;
                if upper {
                    unsafe { tmp = alpha * *a.offset(aj + k) };
                } else {
                    unsafe { tmp = alpha * *a.offset(ak + j) };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *b.offset(bk + i) };
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let ak = k * lda;
                let bk = k * ldb;
                if upper {
                    unsafe { tmp = alpha * *a.offset(ak + j) };
                } else {
                    unsafe { tmp = alpha * *a.offset(aj + k) };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *b.offset(bk + i) };
                    i += 1;
                }
                k += 1;
            }
            j += 1;
        }
    }
}

/// SYR2K  performs one of the symmetric rank 2k operations
/// C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C,
/// where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
/// and  A and B  are  n by k  matrices  in the  first  case  and  k by n matrices in the second case.
/// This is [CSYR2K](http://www.netlib.org/lapack/explore-html/de/d7e/csyr2k_8f.html) and [ZSYR2K](http://www.netlib.org/lapack/explore-html/df/d20/zsyr2k_8f.html) comined in one function
#[inline]
pub fn syr2k<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: isize,
    k: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *const Complex<T>,
    ldb: isize,
    beta: Complex<T>,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
    let upper = uplo == 'u' || uplo == 'U';
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
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldb < max(1, nrowa) {
        info = 9;
    } else if ldc < max(1, n) {
        info = 12;
    }

    if info != 0 {
        panic!("syr2k {}", info);
    }

    let alpha_is_zero = alpha.is_zero();
    let beta_is_one = beta.is_one();
    let beta_is_zero = beta.is_zero();

    if n == 0 || ((alpha_is_zero || k == 0) && beta_is_one) {
        return;
    }

    if alpha_is_zero {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                if beta_is_zero {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                } else {
                    unsafe { *c.offset(cj + i) *= beta };
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            if beta_is_zero {
                let mut i = start;
                while i < stop {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i < stop {
                    unsafe { *c.offset(cj + i) = beta * *c.offset(cj + i) };
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let bl = l * ldb;
                let mut tmp = unsafe { *b.offset(bl + j) };
                let mut tmp2 = unsafe { *a.offset(al + j) };
                if !tmp.is_zero() || !tmp2.is_zero() {
                    tmp *= alpha;
                    tmp2 *= alpha;
                    let mut i = start;
                    while i < stop {
                        unsafe {
                            *c.offset(cj + i) += *a.offset(al + i) * tmp + *b.offset(bl + i) * tmp2
                        };
                        i += 1;
                    }
                }
                l += 1;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let bj = j * ldb;
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                let ai = i * lda;
                let bi = i * ldb;
                let mut tmp = Complex::zero();
                let mut tmp2 = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe {
                        tmp += *a.offset(ai + l) * *b.offset(bj + l);
                        tmp2 += *a.offset(aj + l) * *b.offset(bi + l);
                    }
                    l += 1;
                }
                tmp = alpha * tmp + alpha * tmp2;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
    }
}

/// CSYRK  performs one of the symmetric rank k operations
/// C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C,
/// where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
/// and  A  is an  n by k  matrix in the first case and a  k by n  matrix in the second case.
/// This is [CSYRK](http://www.netlib.org/lapack/explore-html/d3/d6a/csyrk_8f.html) and [ZSYRK](http://www.netlib.org/lapack/explore-html/de/d54/zsyrk_8f.html) comined in one function
#[inline]
pub fn syrk<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: isize,
    k: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    beta: Complex<T>,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let upper = uplo == 'u' || uplo == 'U';
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
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
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldc < max(1, n) {
        info = 10;
    }

    if info != 0 {
        panic!("syrk {}", info);
    }

    let alpha_is_zero = alpha.is_zero();
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();

    if n == 0 || ((alpha_is_zero || k == 0) && beta_is_one) {
        return;
    }

    if alpha_is_zero {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                if beta_is_zero {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                } else {
                    unsafe { *c.offset(cj + i) *= beta };
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            if beta_is_zero {
                let mut i = start;
                while i < stop {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i < stop {
                    unsafe { *c.offset(cj + i) *= beta };
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let mut tmp = unsafe { *a.offset(al + j) };
                if !tmp.is_zero() {
                    tmp *= alpha;
                    let mut i = start;
                    while i < stop {
                        unsafe { *c.offset(cj + i) += tmp * *a.offset(al + i) };
                        i += 1;
                    }
                }
                l += 1;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let aj = j * lda;
            let cj = j * ldc;
            let mut i = start;
            while i < stop {
                let mut tmp = Complex::zero();
                let ai = i * lda;
                let mut l = 0;
                while l < k {
                    unsafe { tmp += *a.offset(ai + l) * *a.offset(aj + l) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    unsafe { tmp += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
    }
}

/// TRMM  performs one of the matrix-matrix operations
///  B := alpha*op( A )*B,   or   B := alpha*B*op( A )
/// where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
/// op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
/// This is [CTRMM](http://www.netlib.org/lapack/explore-html/d4/d9b/ctrmm_8f.html) and [ZTRMM](http://www.netlib.org/lapack/explore-html/d8/de1/ztrmm_8f.html) comined in one function
#[inline]
pub fn trmm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    trans: char,
    diag: char,
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *mut Complex<T>,
    ldb: isize,
) {
    let lside = side == 'l' || side == 'L';

    let nrowa = if lside { m } else { n };
    let nounit = diag == 'n' || diag == 'N';
    let upper = uplo == 'u' || uplo == 'U';
    let alpha_is_zero = alpha.is_zero();
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 3;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 4;
    } else if m < 0 {
        info = 5;
    } else if n < 0 {
        info = 6;
    } else if lda < max(1, nrowa) {
        info = 9;
    } else if ldb < max(1, m) {
        info = 11;
    }
    if info != 0 {
        panic!("trmm {}", info);
    }

    if m == 0 || n == 0 {
        return;
    }

    let noconj = trans == 't' || trans == 'T';

    if alpha_is_zero {
        zero(b, ldb, n, m);
        return;
    }

    if (side == 'l' || side == 'L') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut k = 0;
                while k < m {
                    let mut tmp = unsafe { *b.offset(bj + k) };
                    if !tmp.is_zero() {
                        let ak = k * lda;
                        tmp *= alpha;
                        let mut i = 0;
                        while i < k {
                            unsafe { *b.offset(bj + i) += tmp * *a.offset(ak + i) };
                            i += 1;
                        }
                        if nounit {
                            unsafe { tmp *= *a.offset(ak + k) };
                        }
                        unsafe { *b.offset(bj + k) = tmp };
                    }
                    k += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut k = m;
                while k >= 1 {
                    k -= 1;
                    let mut tmp = unsafe { *b.offset(bj + k) };
                    if !tmp.is_zero() {
                        let ak = k * lda;
                        tmp *= alpha;
                        unsafe { *b.offset(bj + k) = tmp };
                        if nounit {
                            unsafe { *b.offset(bj + k) *= *a.offset(ak + k) };
                        }
                        let mut i = k + 1;
                        while i < m {
                            unsafe { *b.offset(bj + i) += tmp * *a.offset(ak + i) };
                            i += 1;
                        }
                    }
                }
                j += 1;
            }
        }
        return;
    } else if (side == 'l' || side == 'L')
        && (trans == 't' || trans == 'T' || trans == 'c' || trans == 'C')
    {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = m;
                while i >= 1 {
                    i -= 1;
                    let ai = i * lda;
                    let mut tmp = unsafe { *b.offset(bj + i) };
                    if nounit {
                        if noconj {
                            unsafe { tmp *= *a.offset(ai + i) };
                        } else {
                            unsafe { tmp *= a.offset(ai + i).read().conj() };
                        }
                    }
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            unsafe { tmp += *a.offset(ai + k) * *b.offset(bj + k) };
                        } else {
                            unsafe {
                                tmp += a.offset(ai + k).read().conj() * *b.offset(bj + k);
                            };
                        }
                        k += 1;
                    }
                    unsafe { *b.offset(bj + i) = alpha * tmp };
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let ai = i * lda;
                    let mut tmp = unsafe { *b.offset(bj + i) };
                    if nounit {
                        if noconj {
                            unsafe { tmp *= *a.offset(ai + i) };
                        } else {
                            unsafe { tmp *= a.offset(ai + i).read().conj() };
                        }
                    }
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            unsafe { tmp += *a.offset(ai + k) * *b.offset(bj + k) };
                        } else {
                            unsafe { tmp += a.offset(ai + k).read().conj() * *b.offset(bj + k) };
                        }
                        k += 1;
                    }
                    unsafe { *b.offset(bj + i) = alpha * tmp };

                    i += 1;
                }
                j += 1;
            }
        }
        return;
    } else if (side == 'r' || side == 'R') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let aj = j * lda;
                let bj = j * ldb;
                let mut tmp = alpha;
                if nounit {
                    unsafe { tmp *= *a.offset(aj + j) };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bj + i) *= tmp };
                    i += 1;
                }
                let mut k = 0;
                while k < j {
                    let bk = k * ldb;
                    let mut tmp = unsafe { *a.offset(aj + k) };
                    if !tmp.is_zero() {
                        tmp *= alpha;
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) += tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    k += 1;
                }
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let aj = j * lda;
                let mut tmp = alpha;
                if nounit {
                    unsafe { tmp *= *a.offset(aj + j) };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bj + i) *= tmp };
                    i += 1;
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    let tmp2 = unsafe { *a.offset(aj + k) };
                    if !tmp2.is_zero() {
                        tmp = alpha * tmp2;
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) += tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    k += 1;
                }
                j += 1;
            }
        }
        return;
    } else {
        if upper {
            let mut k = 0;
            while k < n {
                let ak = k * lda;
                let bk = k * ldb;
                let mut j = 0;
                while j < k {
                    let bj = j * ldb;
                    let mut tmp = unsafe { *a.offset(ak + j) };
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        }
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) += tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        unsafe { tmp *= *a.offset(ak + k) };
                    } else {
                        unsafe { tmp *= a.offset(ak + k).read().conj() };
                    }
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bk + i) = tmp * *b.offset(bk + i) };
                        i += 1;
                    }
                }
                k += 1;
            }
        } else {
            let mut k = n;
            while k >= 1 {
                k -= 1;
                let ak = k * lda;
                let bk = k * ldb;
                let mut j = k + 1;
                while j < n {
                    let bj = j * lda;
                    let mut tmp = unsafe { *a.offset(ak + j) };
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        };
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) += tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        unsafe { tmp *= *a.offset(ak + k) };
                    } else {
                        unsafe { tmp *= a.offset(ak + k).read().conj() };
                    };
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bk + i) = tmp * *b.offset(bk + i) };
                        i += 1;
                    }
                }
            }
        }
        return;
    }
}

/// TRSM  solves one of the matrix equations
/// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
/// where alpha is a scalar, X and B are m by n matrices, A is a unit, or  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
/// op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
/// The matrix X is overwritten on B.
/// This is [CTRSM](http://www.netlib.org/lapack/explore-html/de/d30/ctrsm_8f.html) and [ZTRSM](http://www.netlib.org/lapack/explore-html/d1/d39/ztrsm_8f.html) comined in one function
#[inline]
pub fn trsm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    trans: char,
    diag: char,
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *mut Complex<T>,
    ldb: isize,
) {
    let lside = side == 'l' || side == 'L';

    let nrowa = if lside { m } else { n };
    let noconj = trans == 't' || trans == 'T';
    let nounit = diag == 'n' || diag == 'N';
    let upper = uplo == 'u' || uplo == 'U';
    let alpha_is_zero = alpha.is_zero();
    let alpha_is_one = alpha.is_one();
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if trans != 'c'
        && trans != 'C'
        && trans != 'n'
        && trans != 'N'
        && trans != 't'
        && trans != 'T'
    {
        info = 3;
    } else if diag != 'n' && diag != 'N' && diag != 'u' && diag != 'U' {
        info = 4;
    } else if m < 0 {
        info = 5;
    } else if n < 0 {
        info = 6;
    } else if lda < max(1, nrowa) {
        info = 9;
    } else if ldb < max(1, m) {
        info = 11;
    }
    if info != 0 {
        panic!("trsm {}", info);
    }

    if m == 0 || n == 0 {
        return;
    }

    if alpha_is_zero {
        zero(b, ldb, n, m);
        return;
    }

    if (side == 'l' || side == 'L') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) = alpha * *b.offset(bj + i) };
                        i += 1;
                    }
                }
                let mut k = m;
                while k >= 1 {
                    k -= 1;
                    let ak = k * lda;
                    let tmp = unsafe { *b.offset(bj + k) };
                    if !tmp.is_zero() {
                        if nounit {
                            unsafe { *b.offset(bj + k) = tmp / a.offset(ak + k).read() };
                        }
                        let tmp = unsafe { *b.offset(bj + k) };
                        let mut i = 0;
                        while i < k {
                            unsafe { *b.offset(bj + i) -= tmp * *a.offset(ak + i) };
                            i += 1;
                        }
                    }
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) = alpha * *b.offset(bj + i) };
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < m {
                    let ak = k * lda;
                    let tmp = unsafe { *b.offset(bj + k) };
                    if !tmp.is_zero() {
                        if nounit {
                            unsafe { *b.offset(bj + k) = tmp / a.offset(ak + k).read() };
                        }
                        let tmp = unsafe { *b.offset(bj + k) };
                        let mut i = k + 1;
                        while i < m {
                            unsafe { *b.offset(bj + i) -= tmp * *a.offset(ak + i) };
                            i += 1;
                        }
                    }
                    k += 1;
                }
                j += 1;
            }
        }
        return;
    }

    if (side == 'l' || side == 'L') && (trans != 'n' && trans != 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let ai = i * lda;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            unsafe { tmp -= *a.offset(ai + k) * *b.offset(bj + k) };
                        } else {
                            unsafe {
                                tmp -= a.offset(ai + k).read().conj() * *b.offset(bj + k);
                            };
                        }
                        k += 1;
                    }
                    if nounit {
                        unsafe { tmp /= *a.offset(ai + i) };
                    }
                    unsafe { *b.offset(bj + i) = tmp };
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = m;
                while i >= 1 {
                    i -= 1;
                    let ai = i * lda;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            unsafe { tmp -= *a.offset(ai + k) * *b.offset(bj + k) };
                        } else {
                            unsafe {
                                tmp -= a.offset(ai + k).read().conj() * *b.offset(bj + k);
                            };
                        }

                        k += 1;
                    }
                    if nounit {
                        if noconj {
                            unsafe { tmp /= *a.offset(ai + i) };
                        } else {
                            unsafe { tmp /= a.offset(ai + i).read().conj() };
                        }
                    }
                    unsafe { *b.offset(bj + i) = tmp };
                }
                j += 1;
            }
        }
        return;
    }

    if (side == 'r' || side == 'R') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let aj = j * lda;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) = alpha * *b.offset(bj + i) };
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < j {
                    let bk = k * ldb;
                    let tmp = unsafe { *a.offset(aj + k) };
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) -= tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = unsafe { Complex::new(T::one(), T::zero()) / *a.offset(aj + j) };
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) *= tmp };
                        i += 1;
                    }
                }
                j += 1;
            }
        } else {
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let bj = j * ldb;
                let aj = j * lda;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) = alpha * *b.offset(bj + i) };
                        i += 1;
                    }
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    let tmp = unsafe { *a.offset(aj + k) };
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            unsafe { *b.offset(bj + i) -= tmp * *b.offset(bk + i) };
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = unsafe { Complex::new(T::one(), T::zero()) / *a.offset(aj + j) };
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) *= tmp };
                        i += 1;
                    }
                }
            }
        }
        return;
    }

    if upper {
        let mut k = n;
        while k >= 1 {
            k -= 1;
            let ak = k * lda;
            let bk = k * ldb;
            if nounit {
                let mut tmp = Complex::new(T::one(), T::zero());
                if noconj {
                    unsafe { tmp /= *a.offset(ak + k) };
                } else {
                    unsafe { tmp /= a.offset(ak + k).read().conj() };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bk + i) = tmp * *b.offset(bk + i) };
                    i += 1;
                }
            }
            let mut j = 0;
            while j < k {
                let bj = j * ldb;
                let mut tmp = unsafe { *a.offset(ak + j) };
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) -= tmp * *b.offset(bk + i) };
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bk + i) = alpha * *b.offset(bk + i) };
                    i += 1;
                }
            }
        }
        return;
    } else {
        let mut k = 0;
        while k < n {
            let ak = k * lda;
            let bk = k * ldb;
            if nounit {
                let mut tmp = Complex::new(T::one(), T::zero());
                if noconj {
                    unsafe { tmp /= *a.offset(ak + k) };
                } else {
                    unsafe { tmp /= a.offset(ak + k).read().conj() };
                }
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bk + i) = tmp * *b.offset(bk + i) };
                    i += 1;
                }
            }
            let mut j = k + 1;
            while j < n {
                let bj = j * ldb;
                let mut tmp = unsafe { *a.offset(ak + j) };
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        unsafe { *b.offset(bj + i) -= tmp * *b.offset(bk + i) };
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    unsafe { *b.offset(bk + i) = alpha * *b.offset(bk + i) };
                    i += 1;
                }
            }
            k += 1;
        }
        return;
    }
}

/// HEMM  performs one of the matrix-matrix operations
/// C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
/// here alpha and beta are scalars, A is an hermitian matrix and  B and C are m by n matrices.
/// This is [CHEMM](http://www.netlib.org/lapack/explore-html/d3/d66/chemm_8f.html) and [ZHEMM](http://www.netlib.org/lapack/explore-html/d6/d3e/zhemm_8f.html) comined in one function
#[inline]
pub fn hemm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    m: isize,
    n: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *const Complex<T>,
    ldb: isize,
    beta: Complex<T>,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let beta_is_zero = beta.is_zero();
    let beta_is_one = beta.is_one();
    let alpha_is_zero = alpha.is_zero();
    let nrowa = if side == 'l' || side == 'L' { m } else { n };
    let upper = uplo == 'u' || uplo == 'U';
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldb < max(1, m) {
        info = 9;
    } else if ldc < max(1, m) {
        info = 12;
    }
    if info != 0 {
        panic!("hemm {}", info);
    }

    if m == 0 || n == 0 || (alpha_is_zero && beta_is_one) {
        return;
    }

    if alpha_is_zero {
        multiply(c, beta, ldc, n, m);
        return;
    }

    if side == 'l' || side == 'L' {
        if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let cj = j * ldc;
                let mut i = 0;
                while i < m {
                    let ai = i * lda;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        unsafe {
                            *c.offset(cj + k) += tmp * *a.offset(ai + k);
                            tmp2 += *b.offset(bj + k) * a.offset(ai + k).read().conj();;
                        }
                        k += 1;
                    }
                    tmp = unsafe { tmp * a.offset(ai + i).read().re + (alpha * tmp2) };
                    if !beta_is_zero {
                        unsafe { tmp += beta * *c.offset(cj + i) };
                    }
                    unsafe { *c.offset(cj + i) = tmp };
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let cj = j * ldc;
                let mut i = m;
                while i >= 1 {
                    i -= 1;
                    let ai = i * lda;
                    let mut tmp = unsafe { alpha * *b.offset(bj + i) };
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        unsafe {
                            *c.offset(cj + k) += tmp * *a.offset(ai + k);
                            tmp2 += *b.offset(bj + k) * a.offset(ai + k).read().conj();
                        }
                        k += 1;
                    }
                    tmp = unsafe { tmp * a.offset(ai + i).read().re + (alpha * tmp2) };
                    if !beta_is_zero {
                        unsafe { tmp += beta * *c.offset(cj + i) };
                    }
                    unsafe { *c.offset(cj + i) = tmp };
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let bj = j * ldb;
            let cj = j * ldc;
            let mut tmp = unsafe { alpha * a.offset(aj + j).read().re };
            let mut i = 0;
            while i < m {
                let mut tmp2 = unsafe { tmp * *b.offset(bj + i) };
                if !beta_is_zero {
                    unsafe { tmp2 += beta * *c.offset(cj + i) };
                }
                unsafe { *c.offset(cj + i) = tmp2 };
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let bk = k * ldb;
                let ak = k * lda;
                if upper {
                    tmp = unsafe { alpha * *a.offset(aj + k) };
                } else {
                    tmp = unsafe { alpha * a.offset(ak + j).read().conj() };
                };
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *b.offset(bk + i) };
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let bk = k * ldb;
                let ak = k * lda;
                if upper {
                    tmp = unsafe { alpha * a.offset(ak + j).read().conj() };
                } else {
                    tmp = unsafe { alpha * *a.offset(aj + k) };
                };
                let mut i = 0;
                while i < m {
                    unsafe { *c.offset(cj + i) += tmp * *b.offset(bk + i) };
                    i += 1;
                }
                k += 1;
            }
            j += 1;
        }
    }
}

/// CHER2K  performs one of the hermitian rank 2k operations
/// C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C, or C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
/// where  alpha and beta  are scalars with  beta  real,  C is an  n by n
/// hermitian matrix and  A and B  are  n by k matrices in the first case and  k by n  matrices in the second case.
/// This is [CHER2K](http://www.netlib.org/lapack/explore-html/d1/d82/cher2k_8f.html) and [ZHER2K](http://www.netlib.org/lapack/explore-html/d7/dfa/zher2k_8f.html) comined in one function
#[inline]
pub fn her2k<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: isize,
    k: isize,
    alpha: Complex<T>,
    a: *const Complex<T>,
    lda: isize,
    b: *const Complex<T>,
    ldb: isize,
    beta: T,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
    let upper = uplo == 'u' || uplo == 'U';
    let alpha_is_zero = alpha.is_zero();
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'n' && trans != 'N' && trans != 'c' && trans != 'C' {
        info = 2;
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldb < max(1, nrowa) {
        info = 9;
    } else if ldc < max(1, n) {
        info = 12;
    }

    if info != 0 {
        panic!("her2k {}", info);
    }

    if n == 0 || ((alpha_is_zero || k == 0) && beta.is_one()) {
        return;
    }

    if alpha_is_zero {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let f1 = if upper { 0 } else { j };
            let f2 = if upper { j + 1 } else { n };
            if beta.is_zero() {
                let mut i = f1;
                while i < f2 {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            }
            let start = if upper { 0 } else { j + 1 };
            let stop = if upper { j - 1 } else { n - 1 };
            let mut i = start;
            while i <= stop {
                unsafe { *c.offset(cj + i) *= beta };
                i += 1;
            }
            let mut ct = unsafe { *c.offset(cj + j) * beta };
            ct.im = T::zero();
            unsafe { *c.offset(cj + j) = ct };
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = 0;
            while j < n {
                let cj = j * ldc;
                if beta.is_zero() {
                    let mut i = 0;
                    while i <= j {
                        unsafe { *c.offset(cj + i) = Complex::zero() };
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < j {
                        unsafe { *c.offset(cj + i) *= beta };
                        i += 1;
                    }
                    let mut ct = unsafe { *c.offset(cj + j) * beta };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + j) = ct };
                } else {
                    let mut ct = unsafe { *c.offset(cj + j) };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + j) = ct };
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    let tmp = unsafe { *b.offset(bl + j) };
                    let tmp2 = unsafe { *a.offset(al + j) };
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        let tmp = alpha * tmp.conj();
                        let tmp2 = (alpha * tmp2).conj();
                        let mut i = 0;
                        while i <= j {
                            unsafe {
                                *c.offset(cj + i) +=
                                    *a.offset(al + i) * tmp + *b.offset(bl + i) * tmp2
                            };
                            i += 1;
                        }
                        let mut ct = unsafe { *c.offset(cj + j) };
                        ct.im = T::zero();
                        unsafe { *c.offset(cj + j) = ct };
                    }
                    l += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let cj = j * ldc;
                if beta.is_zero() {
                    let mut i = j;
                    while i < n {
                        unsafe { *c.offset(cj + i) = Complex::zero() };
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = j + 1;
                    while i < n {
                        unsafe { *c.offset(cj + i) *= beta };
                        i += 1;
                    }
                    let mut ct = unsafe { *c.offset(cj + j) * beta };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + j) = ct };
                } else {
                    let mut ct = unsafe { *c.offset(cj + j) };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + j) = ct };
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    let mut tmp = unsafe { *b.offset(bl + j) };
                    let mut tmp2 = unsafe { *a.offset(al + j) };
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        tmp = alpha * tmp.conj();
                        tmp2 = (alpha * tmp2).conj();
                        let mut i = j + 1;
                        while i < n {
                            unsafe {
                                *c.offset(cj + i) +=
                                    *a.offset(al + i) * tmp + *b.offset(bl + i) * tmp2
                            };
                            i += 1;
                        }
                        let mut ct = unsafe {
                            *c.offset(cj + j)
                                + (*a.offset(al + j) * tmp)
                                + (*b.offset(bl + j) * tmp2)
                        };
                        ct.im = T::zero();
                        unsafe { *c.offset(cj + j) = ct }
                    }
                    l += 1;
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let bj = j * ldb;
            let aj = j * lda;
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                let ai = i * lda;
                let bi = i * ldb;
                let mut tmp = Complex::zero();
                let mut tmp2 = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe {
                        tmp += a.offset(ai + l).read().conj() * *b.offset(bj + l);
                        tmp2 += b.offset(bi + l).read().conj() * *a.offset(aj + l);
                    }
                    l += 1;
                }

                tmp = alpha * tmp + alpha.conj() * tmp2;
                if i == j {
                    tmp.im = T::zero();
                    if beta.is_zero() {
                        unsafe { *c.offset(cj + i) = tmp };
                    } else {
                        unsafe {
                            let mut ct = *c.offset(cj + i);
                            ct.im = T::zero();
                            ct *= beta;
                            *c.offset(cj + i) = ct + tmp;
                        }
                    }
                } else if beta.is_zero() {
                    unsafe { *c.offset(cj + i) = tmp };
                } else {
                    unsafe { *c.offset(cj + i) = *c.offset(cj + i) * beta + tmp };
                }
                i += 1;
            }
            j += 1;
        }
    }
}

/// HERK  performs one of the hermitian rank k operations
/// C := alpha*A*A**H + beta*C, or C := alpha*A**H*A + beta*C,
/// where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
/// matrix and  A  is an  n by k  matrix in the  first case and a  k by n matrix in the second case.
/// This is [CHERK](http://www.netlib.org/lapack/explore-html/d8/d52/cherk_8f.html) and [ZHERK](http://www.netlib.org/lapack/explore-html/d1/db1/zherk_8f.html) comined in one function
#[inline]
pub fn herk<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: isize,
    k: isize,
    alpha: T,
    a: *const Complex<T>,
    lda: isize,
    beta: T,
    c: *mut Complex<T>,
    ldc: isize,
) {
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
    let upper = uplo == 'u' || uplo == 'U';
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'n' && trans != 'N' && trans != 'c' && trans != 'C' {
        info = 2;
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
    } else if lda < max(1, nrowa) {
        info = 7;
    } else if ldc < max(1, n) {
        info = 10;
    }

    if info != 0 {
        panic!("herk {}", info);
    }

    if n == 0 || ((alpha.is_zero() || k == 0) && beta.is_one()) {
        return;
    }

    if alpha.is_zero() {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                if beta.is_zero() {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                } else if j == i {
                    // FIXME remove this line if test passes?
                    let mut ct = unsafe { *c.offset(cj + i) * beta };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + i) = ct };
                } else {
                    unsafe { *c.offset(cj + i) *= beta };
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        let mut j = 0;
        while j < n {
            let fs = if upper { 0 } else { j };
            let fe = if upper { j + 1 } else { n };
            let start = if upper { 0 } else { j + 1 };
            let stop = if upper { j - 1 } else { n - 1 };
            let cj = j * ldc;
            if beta.is_zero() {
                let mut i = fs;
                while i < fe {
                    unsafe { *c.offset(cj + i) = Complex::zero() };
                    i += 1;
                }
            } else if !beta.is_one() {
                let mut i = start;
                while i <= stop {
                    unsafe { *c.offset(cj + i) *= beta };
                    i += 1
                }
                // FIXME only calculate real part.
                let mut ct = unsafe { *c.offset(cj + j) * beta };
                ct.im = T::zero();
                unsafe { *c.offset(cj + j) = ct };
            } else {
                let mut ct = unsafe { *c.offset(cj + j) };
                ct.im = T::zero();
                unsafe { *c.offset(cj + j) = ct };
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let mut tmp = unsafe { *a.offset(al + j) };
                if !tmp.is_zero() {
                    tmp = tmp.conj() * alpha;
                    let mut i = start;
                    while i <= stop {
                        unsafe { *c.offset(cj + i) += tmp * *a.offset(al + i) };
                        i += 1;
                    }
                    let mut ct = unsafe { *c.offset(cj + j) + tmp * *a.offset(al + j) };
                    ct.im = T::zero();
                    unsafe { *c.offset(cj + j) = ct };
                }
                l += 1;
            }
            j += 1;
        }
    } else if upper {
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let cj = j * ldc;
            let mut i = 0;
            while i < j {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe { tmp += a.offset(ai + l).read().conj() * *a.offset(aj + l) };
                    l += 1;
                }
                let mut tmp2 = tmp * alpha;
                if !beta.is_zero() {
                    unsafe { tmp2 += *c.offset(cj + i) * beta };
                }
                unsafe { *c.offset(cj + i) = tmp2 };
                i += 1;
            }
            let mut rtmp = Complex::zero();
            let mut l = 0;
            while l < k {
                unsafe { rtmp += a.offset(aj + l).read().conj() * *a.offset(aj + l) };
                l += 1;
            }
            let mut tmp3: Complex<T> = rtmp * alpha;
            if !beta.is_zero() {
                let re = unsafe { tmp3.re + beta * c.offset(cj + j).read().re };
                tmp3.re = re;
            }
            unsafe { *c.offset(cj + j) = tmp3 };
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let mut rtmp = T::zero();
            let aj = j * lda;
            let cj = j * ldc;
            let mut l = 0;
            while l < k {
                unsafe { rtmp += a.offset(aj + l).read().norm_sqr() };
                l += 1;
            }
            rtmp *= alpha;
            if !beta.is_zero() {
                unsafe { rtmp += beta * c.offset(cj + j).read().re }
            }
            unsafe { *c.offset(cj + j) = Complex::from(rtmp) }
            let mut i = j + 1;
            while i < n {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    unsafe { tmp += a.offset(ai + l).read().conj() * *a.offset(aj + l) };
                    l += 1;
                }
                tmp *= alpha;
                if !beta.is_zero() {
                    unsafe { tmp += *c.offset(cj + i) * beta }
                }
                unsafe { *c.offset(cj + i) = tmp };
                i += 1;
            }
            j += 1;
        }
    }
}
