use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::max;

fn multiply<T: Float + NumAssignOps>(
    left: &mut [Complex<T>],
    right: Complex<T>,
    ld: usize,
    n: usize,
    m: usize,
) {
    if right.is_zero() {
        zero(left, ld, n, m);
    } else {
        let mut j = 0;
        while j < n {
            let coords = j * ld;
            let mut i = 0;
            while i < m {
                left[coords + i] *= right;
                i += 1;
            }
            j += 1;
        }
    }
}

fn zero<T: Float + NumAssignOps>(left: &mut [Complex<T>], ld: usize, n: usize, m: usize) {
    let mut j = 0;
    while j < n {
        let coords = j * ld;
        let mut i = 0;
        while i < m {
            left[coords + i] = Complex::zero();
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
    m: usize,
    n: usize,
    k: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &[Complex<T>],
    ldb: usize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
    ldc: usize,
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
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let tmp = alpha * b[bj + l];
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * a[al + i];
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
    } else if (trans_a == 'n' || trans_a == 'N') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let bl = l * ldb;
                let al = l * lda;
                let tmp = alpha * b[bl + j].conj();
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * a[al + i];
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
    } else if (trans_a == 'n' || trans_a == 'N') && (trans_b == 't' || trans_b == 'T') {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let bl = l * ldb;
                let al = l * lda;
                let tmp = alpha * b[bl + j];
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * a[al + i];
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l].conj() * b[bj + l];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i]
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l].conj() * b[bl + j].conj();
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l].conj() * b[bl + j];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l] * b[bj + l];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l] * b[bl + j].conj();
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
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
                    tmp += a[ai + l] * b[bl + j];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
    }
}

/// CSYMM  performs one of the matrix-matrix operations
/// C := alpha*A*B + beta*C or C := alpha*B*A + beta*C,
/// where  alpha and beta are scalars, A is a symmetric matrix and  B and C are m by n matrices.
/// This is [CSYMM](http://www.netlib.org/lapack/explore-html/db/d59/csymm_8f.html) and [ZSYMM](http://www.netlib.org/lapack/explore-html/df/d51/zsymm_8f.html) comined in one function
#[inline]
pub fn symm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    m: usize,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &[Complex<T>],
    ldb: usize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
    ldc: usize,
) {
    let nrowa = if side == 'l' || side == 'L' { m } else { n };
    let upper = uplo == 'u' || uplo == 'U';
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
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
                    let mut tmp = alpha * b[bj + i];
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k];
                        k += 1
                    }
                    tmp = tmp * a[ai + i] + alpha * tmp2;
                    if !beta_is_zero {
                        tmp += beta * c[cj + i];
                    }
                    c[cj + i] = tmp;
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
                    let mut tmp = alpha * b[bj + i];
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k];
                        k += 1
                    }
                    tmp = tmp * a[ai + i] + alpha * tmp2;
                    if !beta_is_zero {
                        tmp += beta * c[cj + i];
                    }
                    c[cj + i] = tmp;
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
            let mut tmp = alpha * a[aj + j];
            let mut i = 0;
            while i < m {
                let mut tmp2 = tmp * b[bj + i];
                if !beta_is_zero {
                    tmp2 += beta * c[cj + i];
                }
                c[cj + i] = tmp2;
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let ak = k * lda;
                let bk = k * ldb;
                if upper {
                    tmp = alpha * a[aj + k];
                } else {
                    tmp = alpha * a[ak + j];
                }
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * b[bk + i];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let ak = k * lda;
                let bk = k * ldb;
                if upper {
                    tmp = alpha * a[ak + j];
                } else {
                    tmp = alpha * a[aj + k];
                }
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * b[bk + i];
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
    n: usize,
    k: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &[Complex<T>],
    ldb: usize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
    ldc: usize,
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
                    c[cj + i] = Complex::zero();
                } else {
                    c[cj + i] *= beta;
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
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i < stop {
                    c[cj + i] = beta * c[cj + i];
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let bl = l * ldb;
                let mut tmp = b[bl + j];
                let mut tmp2 = a[al + j];
                if !tmp.is_zero() || !tmp2.is_zero() {
                    tmp *= alpha;
                    tmp2 *= alpha;
                    let mut i = start;
                    while i < stop {
                        c[cj + i] += a[al + i] * tmp + b[bl + i] * tmp2;
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
                    tmp += a[ai + l] * b[bj + l];
                    tmp2 += a[aj + l] * b[bi + l];
                    l += 1;
                }
                tmp = alpha * tmp + alpha * tmp2;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
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
    n: usize,
    k: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
    ldc: usize,
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
                    c[cj + i] = Complex::zero();
                } else {
                    c[cj + i] *= beta;
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
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i < stop {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let mut tmp = a[al + j];
                if !tmp.is_zero() {
                    tmp *= alpha;
                    let mut i = start;
                    while i < stop {
                        c[cj + i] += tmp * a[al + i];
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
                    tmp += a[ai + l] * a[aj + l];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[cj + i];
                }
                c[cj + i] = tmp;
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
    m: usize,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &mut [Complex<T>],
    ldb: usize,
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
                    let mut tmp = b[bj + k];
                    if !tmp.is_zero() {
                        let ak = k * lda;
                        tmp *= alpha;
                        let mut i = 0;
                        while i < k {
                            b[bj + i] += tmp * a[ak + i];
                            i += 1;
                        }
                        if nounit {
                            tmp *= a[ak + k];
                        }
                        b[bj + k] = tmp;
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
                    let mut tmp = b[bj + k];
                    if !tmp.is_zero() {
                        let ak = k * lda;
                        tmp *= alpha;
                        b[bj + k] = tmp;
                        if nounit {
                            b[bj + k] *= a[ak + k];
                        }
                        let mut i = k + 1;
                        while i < m {
                            b[bj + i] += tmp * a[ak + i];
                            i += 1;
                        }
                    }
                }
                j += 1;
            }
        }
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
                    let mut tmp = b[bj + i];
                    if nounit {
                        if noconj {
                            tmp *= a[ai + i];
                        } else {
                            tmp *= a[ai + i].conj();
                        }
                    }
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            tmp += a[ai + k] * b[bj + k];
                        } else {
                            tmp += a[ai + k].conj() * b[bj + k];
                        }
                        k += 1;
                    }
                    b[bj + i] = alpha * tmp;
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
                    let mut tmp = b[bj + i];
                    if nounit {
                        if noconj {
                            tmp *= a[ai + i];
                        } else {
                            tmp *= a[ai + i].conj();
                        }
                    }
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            tmp += a[ai + k] * b[bj + k];
                        } else {
                            tmp += a[ai + k].conj() * b[bj + k];
                        }
                        k += 1;
                    }
                    b[bj + i] = alpha * tmp;

                    i += 1;
                }
                j += 1;
            }
        }
    } else if (side == 'r' || side == 'R') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let aj = j * lda;
                let bj = j * ldb;
                let mut tmp = alpha;
                if nounit {
                    tmp *= a[aj + j];
                }
                let mut i = 0;
                while i < m {
                    b[bj + i] *= tmp;
                    i += 1;
                }
                let mut k = 0;
                while k < j {
                    let bk = k * ldb;
                    let mut tmp = a[aj + k];
                    if !tmp.is_zero() {
                        tmp *= alpha;
                        let mut i = 0;
                        while i < m {
                            b[bj + i] += tmp * b[bk + i];
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
                    tmp *= a[aj + j];
                }
                let mut i = 0;
                while i < m {
                    b[bj + i] *= tmp;
                    i += 1;
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    let tmp2 = a[aj + k];
                    if !tmp2.is_zero() {
                        tmp = alpha * tmp2;
                        let mut i = 0;
                        while i < m {
                            b[bj + i] += tmp * b[bk + i];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                j += 1;
            }
        }
    } else {
        if upper {
            let mut k = 0;
            while k < n {
                let ak = k * lda;
                let bk = k * ldb;
                let mut j = 0;
                while j < k {
                    let bj = j * ldb;
                    let mut tmp = a[ak + j];
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        }
                        let mut i = 0;
                        while i < m {
                            b[bj + i] += tmp * b[bk + i];
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        tmp *= a[ak + k];
                    } else {
                        tmp *= a[ak + k].conj();
                    }
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[bk + i] = tmp * b[bk + i];
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
                    let mut tmp = a[ak + j];
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        };
                        let mut i = 0;
                        while i < m {
                            b[bj + i] += tmp * b[bk + i];
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        tmp *= a[ak + k];
                    } else {
                        tmp *= a[ak + k].conj();
                    };
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[bk + i] = tmp * b[bk + i];
                        i += 1;
                    }
                }
            }
        }
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
    m: usize,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &mut [Complex<T>],
    ldb: usize,
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
                        b[bj + i] = alpha * b[bj + i];
                        i += 1;
                    }
                }
                let mut k = m;
                while k >= 1 {
                    k -= 1;
                    let ak = k * lda;
                    let tmp = b[bj + k];
                    if !tmp.is_zero() {
                        if nounit {
                            b[bj + k] = tmp / a[ak + k];
                        }
                        let tmp = b[bj + k];
                        let mut i = 0;
                        while i < k {
                            b[bj + i] -= tmp * a[ak + i];
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
                        b[bj + i] = alpha * b[bj + i];
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < m {
                    let ak = k * lda;
                    let tmp = b[bj + k];
                    if !tmp.is_zero() {
                        if nounit {
                            b[bj + k] = tmp / a[ak + k];
                        }
                        let tmp = b[bj + k];
                        let mut i = k + 1;
                        while i < m {
                            b[bj + i] -= tmp * a[ak + i];
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
                    let mut tmp = alpha * b[bj + i];
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            tmp -= a[ai + k] * b[bj + k];
                        } else {
                            tmp -= a[ai + k].conj() * b[bj + k];
                        }
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[ai + i];
                    }
                    b[bj + i] = tmp;
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
                    let mut tmp = alpha * b[bj + i];
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            tmp -= a[ai + k] * b[bj + k];
                        } else {
                            tmp -= a[ai + k].conj() * b[bj + k];
                        }

                        k += 1;
                    }
                    if nounit {
                        if noconj {
                            tmp /= a[ai + i];
                        } else {
                            tmp /= a[ai + i].conj();
                        }
                    }
                    b[bj + i] = tmp;
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
                        b[bj + i] = alpha * b[bj + i];
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < j {
                    let bk = k * ldb;
                    let tmp = a[aj + k];
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            b[bj + i] -= tmp * b[bk + i];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = Complex::new(T::one(), T::zero()) / a[aj + j];
                    let mut i = 0;
                    while i < m {
                        b[bj + i] *= tmp;
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
                        b[bj + i] = alpha * b[bj + i];
                        i += 1;
                    }
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    let tmp = a[aj + k];
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            b[bj + i] -= tmp * b[bk + i];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = Complex::new(T::one(), T::zero()) / a[aj + j];
                    let mut i = 0;
                    while i < m {
                        b[bj + i] *= tmp;
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
                    tmp /= a[ak + k];
                } else {
                    tmp /= a[ak + k].conj();
                }
                let mut i = 0;
                while i < m {
                    b[bk + i] = tmp * b[bk + i];
                    i += 1;
                }
            }
            let mut j = 0;
            while j < k {
                let bj = j * ldb;
                let mut tmp = a[ak + j];
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        b[bj + i] -= tmp * b[bk + i];
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    b[bk + i] = alpha * b[bk + i];
                    i += 1;
                }
            }
        }
    } else {
        let mut k = 0;
        while k < n {
            let ak = k * lda;
            let bk = k * ldb;
            if nounit {
                let mut tmp = Complex::new(T::one(), T::zero());
                if noconj {
                    tmp /= a[ak + k];
                } else {
                    tmp /= a[ak + k].conj();
                }
                let mut i = 0;
                while i < m {
                    b[bk + i] = tmp * b[bk + i];
                    i += 1;
                }
            }
            let mut j = k + 1;
            while j < n {
                let bj = j * ldb;
                let mut tmp = a[ak + j];
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        b[bj + i] -= tmp * b[bk + i];

                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    b[bk + i] = alpha * b[bk + i];

                    i += 1;
                }
            }
            k += 1;
        }
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
    m: usize,
    n: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &[Complex<T>],
    ldb: usize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
    ldc: usize,
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
        multiply(c, beta, ldc as usize, n as usize, m as usize);
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
                    let mut tmp = alpha * b[bj + i];
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k].conj();
                        k += 1;
                    }
                    tmp = tmp * a[ai + i].re + (alpha * tmp2);
                    if !beta_is_zero {
                        tmp += beta * c[cj + i];
                    }
                    c[cj + i] = tmp;
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
                    let mut tmp = alpha * b[bj + i];
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k].conj();
                        k += 1;
                    }
                    tmp = tmp * a[ai + i].re + (alpha * tmp2);
                    if !beta_is_zero {
                        tmp += beta * c[cj + i];
                    }
                    c[cj + i] = tmp;
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
            let mut tmp = alpha * a[aj + j].re;
            let mut i = 0;
            while i < m {
                let mut tmp2 = tmp * b[bj + i];
                if !beta_is_zero {
                    tmp2 += beta * c[cj + i];
                }
                c[cj + i] = tmp2;
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let bk = k * ldb;
                let ak = k * lda;
                if upper {
                    tmp = alpha * a[aj + k];
                } else {
                    tmp = alpha * a[ak + j].conj();
                };
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * b[bk + i];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let bk = k * ldb;
                let ak = k * lda;
                if upper {
                    tmp = alpha * a[ak + j].conj();
                } else {
                    tmp = alpha * a[aj + k];
                };
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * b[bk + i];
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
    n: usize,
    k: usize,
    alpha: Complex<T>,
    a: &[Complex<T>],
    lda: usize,
    b: &[Complex<T>],
    ldb: usize,
    beta: T,
    c: &mut [Complex<T>],
    ldc: usize,
) {
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
    let upper = uplo == 'u' || uplo == 'U';
    let alpha_is_zero = alpha.is_zero();
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'n' && trans != 'N' && trans != 'c' && trans != 'C' {
        info = 2;
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
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            }
            let start = if upper { 0 as isize } else { j as isize + 1 };
            let stop = if upper {
                j as isize - 1
            } else {
                n as isize - 1
            };
            let mut i = start;
            while i <= stop {
                c[(cj as isize + i) as usize] *= beta;
                i += 1;
            }
            let mut ct = c[cj + j] * beta;
            ct.im = T::zero();
            c[cj + j] = ct;
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
                        c[cj + i] = Complex::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < j {
                        c[cj + i] *= beta;
                        i += 1;
                    }
                    let mut ct = c[cj + j] * beta;
                    ct.im = T::zero();
                    c[cj + j] = ct;
                } else {
                    let mut ct = c[cj + j];
                    ct.im = T::zero();
                    c[cj + j] = ct;
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    let tmp = b[bl + j];
                    let tmp2 = a[al + j];
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        let tmp = alpha * tmp.conj();
                        let tmp2 = (alpha * tmp2).conj();
                        let mut i = 0;
                        while i <= j {
                            c[cj + i] += a[al + i] * tmp + b[bl + i] * tmp2;
                            i += 1;
                        }
                        let mut ct = c[cj + j];
                        ct.im = T::zero();
                        c[cj + j] = ct;
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
                        c[cj + i] = Complex::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = j + 1;
                    while i < n {
                        c[cj + i] *= beta;
                        i += 1;
                    }
                    let mut ct = c[cj + j] * beta;
                    ct.im = T::zero();
                    c[cj + j] = ct;
                } else {
                    let mut ct = c[cj + j];
                    ct.im = T::zero();
                    c[cj + j] = ct;
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    let mut tmp = b[bl + j];
                    let mut tmp2 = a[al + j];
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        tmp = alpha * tmp.conj();
                        tmp2 = (alpha * tmp2).conj();
                        let mut i = j + 1;
                        while i < n {
                            c[cj + i] += a[al + i] * tmp + b[bl + i] * tmp2;
                            i += 1;
                        }
                        let mut ct = c[cj + j] + (a[al + j] * tmp) + (b[bl + j] * tmp2);
                        ct.im = T::zero();
                        c[cj + j] = ct
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
                    tmp += a[ai + l].conj() * b[bj + l];
                    tmp2 += b[bi + l].conj() * a[aj + l];
                    l += 1;
                }

                tmp = alpha * tmp + alpha.conj() * tmp2;
                if i == j {
                    tmp.im = T::zero();
                    if beta.is_zero() {
                        c[cj + i] = tmp;
                    } else {
                        let mut ct = c[cj + i];
                        ct.im = T::zero();
                        ct *= beta;
                        c[cj + i] = ct + tmp;
                    }
                } else if beta.is_zero() {
                    c[cj + i] = tmp;
                } else {
                    c[cj + i] = c[cj + i] * beta + tmp;
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
    n: usize,
    k: usize,
    alpha: T,
    a: &[Complex<T>],
    lda: usize,
    beta: T,
    c: &mut [Complex<T>],
    ldc: usize,
) {
    let nrowa = if trans == 'n' || trans == 'N' { n } else { k };
    let upper = uplo == 'u' || uplo == 'U';
    let mut info = 0;
    if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 1;
    } else if trans != 'n' && trans != 'N' && trans != 'c' && trans != 'C' {
        info = 2;
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
                    c[cj + i] = Complex::zero();
                } else if j == i {
                    // FIXME remove this line if test passes?
                    let mut ct = c[cj + i] * beta;
                    ct.im = T::zero();
                    c[cj + i] = ct;
                } else {
                    c[cj + i] *= beta;
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
            let start = if upper { 0 } else { j as isize + 1 };
            let stop = if upper {
                j as isize - 1
            } else {
                n as isize - 1
            };
            let cj = j * ldc;
            if beta.is_zero() {
                let mut i = fs;
                while i < fe {
                    c[cj + i] = Complex::zero();
                    i += 1;
                }
            } else if !beta.is_one() {
                let mut i = start;
                while i <= stop {
                    c[(cj as isize + i) as usize] *= beta;
                    i += 1
                }
                // FIXME only calculate real part.
                let mut ct = c[cj + j] * beta;
                ct.im = T::zero();
                c[cj + j] = ct;
            } else {
                let mut ct = c[cj + j];
                ct.im = T::zero();
                c[cj + j] = ct;
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                let mut tmp = a[al + j];
                if !tmp.is_zero() {
                    tmp = tmp.conj() * alpha;
                    let mut i = start;
                    while i <= stop {
                        c[(cj as isize + i) as usize] += tmp * a[(al as isize + i) as usize];
                        i += 1;
                    }
                    let mut ct = c[cj + j] + tmp * a[al + j];
                    ct.im = T::zero();
                    c[cj + j] = ct;
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
                    tmp += a[ai + l].conj() * a[aj + l];
                    l += 1;
                }
                let mut tmp2 = tmp * alpha;
                if !beta.is_zero() {
                    tmp2 += c[cj + i] * beta;
                }
                c[cj + i] = tmp2;
                i += 1;
            }
            let mut rtmp = Complex::zero();
            let mut l = 0;
            while l < k {
                rtmp += a[aj + l].conj() * a[aj + l];
                l += 1;
            }
            let mut tmp3: Complex<T> = rtmp * alpha;
            if !beta.is_zero() {
                let re = tmp3.re + beta * c[cj + j].re;
                tmp3.re = re;
            }
            c[cj + j] = tmp3;
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
                rtmp += a[aj + l].norm_sqr();
                l += 1;
            }
            rtmp *= alpha;
            if !beta.is_zero() {
                rtmp += beta * c[cj + j].re
            }
            c[cj + j].re = rtmp;
            c[cj + j].im = T::zero();
            let mut i = j + 1;
            while i < n {
                let ai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[ai + l].conj() * a[aj + l];
                    l += 1;
                }
                tmp *= alpha;
                if !beta.is_zero() {
                    tmp += c[cj + i] * beta
                }
                c[cj + i] = tmp;
                i += 1;
            }
            j += 1;
        }
    }
}
