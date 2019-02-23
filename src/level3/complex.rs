use num_complex::Complex;
use num_traits::{Float, NumAssignOps, One, Zero};
use std::cmp::max;

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
    a: &[Complex<T>],
    lda: isize,
    b: &[Complex<T>],
    ldb: isize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
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
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                if beta_is_zero {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else {
                    c[(coorcj + i) as usize] *= beta;
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if (trans_a == 'n' || trans_a == 'N') && (trans_b == 'n' || trans_b == 'N') {
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let coorbj = j * ldb;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let cooral = l * lda;
                let tmp = alpha * b[(coorbj + l) as usize];
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] += tmp * a[(cooral + i) as usize];
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
            let coorcj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coorbl = l * ldb;
                let cooral = l * lda;
                let tmp = alpha * b[(coorbl + j) as usize].conj();
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] += tmp * a[(cooral + i) as usize];
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
            let coorcj = j * ldc;
            if beta_is_zero {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coorbl = l * ldb;
                let cooral = l * lda;
                let tmp = alpha * b[(coorbl + j) as usize];
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] += tmp * a[(cooral + i) as usize];
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
            let coorbj = j * ldb;
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coorai + l) as usize].conj() * b[(coorbj + l) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize]
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'c' || trans_a == 'C') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let coorbl = l * ldb;
                    tmp += a[(coorai + l) as usize].conj() * b[(coorbl + j) as usize].conj();
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 'c' || trans_a == 'C') && (trans_b == 't' || trans_b == 'T') {
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let coorbl = l * ldb;
                    tmp += a[(coorai + l) as usize].conj() * b[(coorbl + j) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 't' || trans_a == 'T') && (trans_b == 'n' || trans_b == 'N') {
        let mut j = 0;
        while j < n {
            let coorbj = j * ldb;
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coorai + l) as usize] * b[(coorbj + l) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
        return;
    } else if (trans_a == 't' || trans_a == 'T') && (trans_b == 'c' || trans_b == 'C') {
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let coorbl = l * ldb;
                    tmp += a[(coorai + l) as usize] * b[(coorbl + j) as usize].conj();
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
        return;
    } else {
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    let coorbl = l * ldb;
                    tmp += a[(coorai + l) as usize] * b[(coorbl + j) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp;
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
    a: &[Complex<T>],
    lda: isize,
    b: &[Complex<T>],
    ldb: isize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
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
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                if beta_is_zero {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else {
                    c[(coorcj + i) as usize] *= beta;
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if side == 'l' || side == 'L' {
        if upper {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let coord_cj = j * ldc;
                let mut i = 0;
                while i < m {
                    let coord_ai = i * ldc;
                    let mut tmp = alpha * b[(coord_bj + i) as usize];
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        c[(coord_cj + k) as usize] += tmp * a[(coord_ai + k) as usize];
                        tmp2 += b[(coord_bj + k) as usize] * a[(coord_ai + k) as usize];
                        k += 1
                    }
                    tmp = tmp * a[(coord_ai + i) as usize] + alpha * tmp2;
                    if !beta_is_zero {
                        tmp += beta * c[(coord_cj + i) as usize];
                    }
                    c[(coord_cj + i) as usize] = tmp;
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let coord_cj = j * ldc;
                let mut i = m - 1;
                while i >= 0 {
                    let coord_ai = i * lda;
                    let mut tmp = alpha * b[(coord_bj + i) as usize];
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        c[(coord_cj + k) as usize] += tmp * a[(coord_ai + k) as usize];
                        tmp2 += b[(coord_bj + k) as usize] * a[(coord_ai + k) as usize];
                        k += 1
                    }
                    tmp = tmp * a[(coord_ai + i) as usize] + alpha * tmp2;
                    if !beta_is_zero {
                        tmp += beta * c[(coord_cj + i) as usize];
                    }
                    c[(coord_cj + i) as usize] = tmp;
                    i -= 1
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let coord_aj = j * lda;
            let coord_bj = j * ldb;
            let coord_cj = j * ldc;
            let mut tmp = alpha * a[(coord_aj + j) as usize];
            let mut i = 0;
            while i < m {
                let mut tmp2 = tmp * b[(coord_bj + i) as usize];
                if !beta_is_zero {
                    tmp2 += beta * c[(coord_cj + i) as usize];
                }
                c[(coord_cj + i) as usize] = tmp2;
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let coor_ak = k * lda;
                let coor_bk = k * ldb;
                if upper {
                    tmp = alpha * a[(coord_aj + k) as usize];
                } else {
                    tmp = alpha * a[(coor_ak + j) as usize];
                }
                let mut i = 0;
                while i < m {
                    c[(coord_cj + i) as usize] += tmp * b[(coor_bk + i) as usize];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let coor_ak = k * lda;
                let coor_bk = k * ldb;
                if upper {
                    tmp = alpha * a[(coor_ak + j) as usize];
                } else {
                    tmp = alpha * a[(coord_aj + k) as usize];
                }
                let mut i = 0;
                while i < m {
                    c[(coord_cj + i) as usize] += tmp * b[(coor_bk + i) as usize];
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
    a: &[Complex<T>],
    lda: isize,
    b: &[Complex<T>],
    ldb: isize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
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
            let coorcj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                if beta_is_zero {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else {
                    c[(coorcj + i) as usize] *= beta;
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
            let coord_cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            if beta_is_zero {
                let mut i = start;
                while i <= stop {
                    c[(coord_cj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i <= stop {
                    c[(coord_cj + i) as usize] = beta * c[(coord_cj + i) as usize];
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coord_al = l * lda;
                let coord_bl = l * ldb;
                let mut tmp = b[(coord_bl + j) as usize];
                let mut tmp2 = a[(coord_al + j) as usize];
                if !tmp.is_zero() || !tmp2.is_zero() {
                    tmp *= alpha;
                    tmp2 *= alpha;
                    let mut i = start;
                    while i <= stop {
                        c[(coord_cj + i) as usize] +=
                            a[(coord_al + i) as usize] * tmp + b[(coord_bl + i) as usize] * tmp2;
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
            let coord_aj = j * lda;
            let coord_bj = j * ldb;
            let coord_cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                let coor_ai = i * lda;
                let coor_bi = i * ldb;
                let mut tmp = Complex::zero();
                let mut tmp2 = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coor_ai + l) as usize] * b[(coord_bj + l) as usize];
                    tmp2 += a[(coord_aj + l) as usize] * b[(coor_bi + l) as usize];
                    l += 1;
                }
                tmp = alpha * tmp + alpha * tmp2;
                if !beta_is_zero {
                    tmp += beta * c[(coord_cj + i) as usize];
                }
                c[(coord_cj + i) as usize] = tmp;
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
    a: &[Complex<T>],
    lda: isize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
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
            let coorcj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                if beta_is_zero {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else {
                    c[(coorcj + i) as usize] *= beta;
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
            let coord_cj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            if beta_is_zero {
                let mut i = start;
                while i <= stop {
                    c[(coord_cj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta_is_one {
                let mut i = start;
                while i <= stop {
                    c[(coord_cj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coord_al = l * lda;
                let mut tmp = a[(coord_al + j) as usize];
                if !tmp.is_zero() {
                    tmp *= alpha;
                    let mut i = start;
                    while i <= stop {
                        c[(coord_cj + i) as usize] += tmp * a[(coord_al + i) as usize];
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
            let stop = if upper { j } else { n - 1 };
            let coord_aj = j * lda;
            let coord_cj = j * ldc;
            let mut i = start;
            while i <= stop {
                let mut tmp = Complex::zero();
                let coord_ai = i * lda;
                let mut l = 0;
                while l < k {
                    tmp += a[(coord_ai + l) as usize] * a[(coord_aj + l) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta_is_zero {
                    tmp += beta * c[(coord_cj + i) as usize];
                }
                c[(coord_cj + i) as usize] = tmp;
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
    a: &[Complex<T>],
    lda: isize,
    b: &mut [Complex<T>],
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
        let mut j = 0;
        while j < n {
            let coorbj = j * ldb;
            let mut i = 0;
            while i < m {
                b[(coorbj + i) as usize] = Complex::zero();
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if (side == 'l' || side == 'L') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let mut k = 0;
                while k < m {
                    let mut tmp = b[(coorbj + k) as usize];
                    if !tmp.is_zero() {
                        let coorak = k * lda;
                        tmp *= alpha;
                        let mut i = 0;
                        while i < k {
                            b[(coorbj + i) as usize] += tmp * a[(coorak + i) as usize];
                            i += 1;
                        }
                        if nounit {
                            tmp *= a[(coorak + k) as usize];
                        }
                        b[(coorbj + k) as usize] = tmp;
                    }
                    k += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let mut k = m - 1;
                while k >= 0 {
                    let mut tmp = b[(coorbj + k) as usize];
                    if !tmp.is_zero() {
                        let coorak = k * lda;
                        tmp *= alpha;
                        b[(coorbj + k) as usize] = tmp;
                        if nounit {
                            b[(coorbj + k) as usize] *= a[(coorak + k) as usize];
                        }
                        let mut i = k + 1;
                        while i < m {
                            b[(coorbj + i) as usize] += tmp * a[(coorak + i) as usize];
                            i += 1;
                        }
                    }
                    k -= 1;
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
                let coorbj = j * ldb;
                let mut i = m - 1;
                while i >= 0 {
                    let coorai = i * lda;
                    let mut tmp = b[(coorbj + i) as usize];
                    if nounit {
                        if noconj {
                            tmp *= a[(coorai + i) as usize];
                        } else {
                            tmp *= a[(coorai + i) as usize].conj();
                        }
                    }
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            tmp += a[(coorai + k) as usize] * b[(coorbj + k) as usize];
                        } else {
                            tmp += a[(coorai + k) as usize].conj() * b[(coorbj + k) as usize];
                        }
                        k += 1;
                    }
                    b[(coorbj + i) as usize] = alpha * tmp;
                    i -= 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let mut i = 0;
                while i < m {
                    let coorai = i * lda;
                    let mut tmp = b[(coorbj + i) as usize];
                    if nounit {
                        if noconj {
                            tmp *= a[(coorai + i) as usize];
                        } else {
                            tmp *= a[(coorai + i) as usize].conj();
                        }
                    }
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            tmp += a[(coorai + k) as usize] * b[(coorbj + k) as usize];
                        } else {
                            tmp += a[(coorai + k) as usize].conj() * b[(coorbj + k) as usize];
                        }
                        k += 1;
                    }
                    b[(coorbj + i) as usize] = alpha * tmp;

                    i += 1;
                }
                j += 1;
            }
        }
        return;
    } else if (side == 'r' || side == 'R') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = n - 1;
            while j >= 0 {
                let cooraj = j * lda;
                let coorbj = j * ldb;
                let mut tmp = alpha;
                if nounit {
                    tmp *= a[(cooraj + j) as usize];
                }
                let mut i = 0;
                while i < m {
                    b[(coorbj + i) as usize] *= tmp;
                    i += 1;
                }
                let mut k = 0;
                while k < j {
                    let coorbk = k * ldb;
                    let mut tmp = a[(cooraj + k) as usize];
                    if !tmp.is_zero() {
                        tmp *= alpha;
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] += tmp * b[(coorbk + i) as usize];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                j -= 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let cooraj = j * lda;
                let mut tmp = alpha;
                if nounit {
                    tmp *= a[(cooraj + j) as usize];
                }
                let mut i = 0;
                while i < m {
                    b[(coorbj + i) as usize] *= tmp;
                    i += 1;
                }
                let mut k = j + 1;
                while k < n {
                    let coorbk = k * ldb;
                    let tmp2 = a[(cooraj + k) as usize];
                    if !tmp2.is_zero() {
                        tmp = alpha * tmp2;
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] += tmp * b[(coorbk + i) as usize];
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
                let coorak = k * lda;
                let coorbk = k * ldb;
                let mut j = 0;
                while j < k {
                    let coorbj = j * ldb;
                    let mut tmp = a[(coorak + j) as usize];
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        }
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] += tmp * b[(coorbk + i) as usize];
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        tmp *= a[(coorak + k) as usize];
                    } else {
                        tmp *= a[(coorak + k) as usize].conj();
                    }
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[(coorbk + i) as usize] = tmp * b[(coorbk + i) as usize];
                        i += 1;
                    }
                }
                k += 1;
            }
        } else {
            let mut k = n - 1;
            while k >= 0 {
                let coorak = k * lda;
                let coorbk = k * ldb;
                let mut j = k + 1;
                while j < n {
                    let coorbj = j * lda;
                    let mut tmp = a[(coorak + j) as usize];
                    if !tmp.is_zero() {
                        if noconj {
                            tmp *= alpha;
                        } else {
                            tmp = alpha * tmp.conj();
                        };
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] += tmp * b[(coorbk + i) as usize];
                            i += 1;
                        }
                    }
                    j += 1;
                }
                let mut tmp = alpha;
                if nounit {
                    if noconj {
                        tmp *= a[(coorak + k) as usize];
                    } else {
                        tmp *= a[(coorak + k) as usize].conj();
                    };
                }
                if !tmp.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[(coorbk + i) as usize] = tmp * b[(coorbk + i) as usize];
                        i += 1;
                    }
                }
                k -= 1;
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
    a: &[Complex<T>],
    lda: isize,
    b: &mut [Complex<T>],
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
        let mut j = 0;
        while j < n {
            let coorbj = j * ldb;
            let mut i = 0;
            while i < m {
                b[(coorbj + i) as usize] = Complex::zero();
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if (side == 'l' || side == 'L') && (trans == 'n' || trans == 'N') {
        if upper {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] = alpha * b[(coorbj + i) as usize];
                        i += 1;
                    }
                }
                let mut k = m - 1;
                while k >= 0 {
                    let coorak = k * lda;
                    let tmp = b[(coorbj + k) as usize];
                    if !tmp.is_zero() {
                        if nounit {
                            b[(coorbj + k) as usize] = tmp / a[(coorak + k) as usize];
                        }
                        let tmp = b[(coorbj + k) as usize];
                        let mut i = 0;
                        while i < k {
                            b[(coorbj + i) as usize] -= tmp * a[(coorak + i) as usize];
                            i += 1;
                        }
                    }
                    k -= 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] = alpha * b[(coorbj + i) as usize];
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < m {
                    let coorak = k * lda;
                    let tmp = b[(coorbj + k) as usize];
                    if !tmp.is_zero() {
                        if nounit {
                            b[(coorbj + k) as usize] = tmp / a[(coorak + k) as usize];
                        }
                        let tmp = b[(coorbj + k) as usize];
                        let mut i = k + 1;
                        while i < m {
                            b[(coorbj + i) as usize] -= tmp * a[(coorak + i) as usize];
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
                let coorbj = j * ldb;
                let mut i = 0;
                while i < m {
                    let coorai = i * lda;
                    let mut tmp = alpha * b[(coorbj + i) as usize];
                    let mut k = 0;
                    while k < i {
                        if noconj {
                            tmp -= a[(coorai + k) as usize] * b[(coorbj + k) as usize];
                        } else {
                            tmp -= a[(coorai + k) as usize].conj() * b[(coorbj + k) as usize];
                        }
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[(coorai + i) as usize];
                    }
                    b[(coorbj + i) as usize] = tmp;
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let mut i = m - 1;
                while i >= 0 {
                    let coorai = i * lda;
                    let mut tmp = alpha * b[(coorbj + i) as usize];
                    let mut k = i + 1;
                    while k < m {
                        if noconj {
                            tmp -= a[(coorai + k) as usize] * b[(coorbj + k) as usize];
                        } else {
                            tmp -= a[(coorai + k) as usize].conj() * b[(coorbj + k) as usize];
                        }

                        k += 1;
                    }
                    if nounit {
                        if noconj {
                            tmp /= a[(coorai + i) as usize];
                        } else {
                            tmp /= a[(coorai + i) as usize].conj();
                        }
                    }
                    b[(coorbj + i) as usize] = tmp;
                    i -= 1;
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
                let coorbj = j * ldb;
                let cooraj = j * lda;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] = alpha * b[(coorbj + i) as usize];
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < j {
                    let coorbk = k * ldb;
                    let tmp = a[(cooraj + k) as usize];
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] -= tmp * b[(coorbk + i) as usize];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = Complex::new(T::one(), T::zero()) / a[(cooraj + j) as usize];
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] *= tmp;
                        i += 1;
                    }
                }
                j += 1;
            }
        } else {
            let mut j = n - 1;
            while j >= 0 {
                let coorbj = j * ldb;
                let cooraj = j * lda;
                if !alpha_is_one {
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] = alpha * b[(coorbj + i) as usize];
                        i += 1;
                    }
                }
                let mut k = j + 1;
                while k < n {
                    let coorbk = k * ldb;
                    let tmp = a[(cooraj + k) as usize];
                    if !tmp.is_zero() {
                        let mut i = 0;
                        while i < m {
                            b[(coorbj + i) as usize] -= tmp * b[(coorbk + i) as usize];
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = Complex::new(T::one(), T::zero()) / a[(cooraj + j) as usize];
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] *= tmp;
                        i += 1;
                    }
                }
                j -= 1;
            }
        }
        return;
    }

    if upper {
        let mut k = n - 1;
        while k >= 0 {
            let coorak = k * lda;
            let coorbk = k * ldb;
            if nounit {
                let mut tmp = Complex::new(T::one(), T::zero());
                if noconj {
                    tmp /= a[(coorak + k) as usize];
                } else {
                    tmp /= a[(coorak + k) as usize].conj();
                }
                let mut i = 0;
                while i < m {
                    b[(coorbk + i) as usize] = tmp * b[(coorbk + i) as usize];
                    i += 1;
                }
            }
            let mut j = 0;
            while j < k {
                let coorbj = j * ldb;
                let mut tmp = a[(coorak + j) as usize];
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] -= tmp * b[(coorbk + i) as usize];
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    b[(coorbk + i) as usize] = alpha * b[(coorbk + i) as usize];
                    i += 1;
                }
            }
            k -= 1;
        }
        return;
    } else {
        let mut k = 0;
        while k < n {
            let coorak = k * lda;
            let coorbk = k * ldb;
            if nounit {
                let mut tmp = Complex::new(T::one(), T::zero());
                if noconj {
                    tmp /= a[(coorak + k) as usize];
                } else {
                    tmp /= a[(coorak + k) as usize].conj();
                }
                let mut i = 0;
                while i < m {
                    b[(coorbk + i) as usize] = tmp * b[(coorbk + i) as usize];
                    i += 1;
                }
            }
            let mut j = k + 1;
            while j < n {
                let coorbj = j * ldb;
                let mut tmp = a[(coorak + j) as usize];
                if !tmp.is_zero() {
                    if !noconj {
                        tmp = tmp.conj()
                    }
                    let mut i = 0;
                    while i < m {
                        b[(coorbj + i) as usize] -= tmp * b[(coorbk + i) as usize];

                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha_is_one {
                let mut i = 0;
                while i < m {
                    b[(coorbk + i) as usize] = alpha * b[(coorbk + i) as usize];

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
    a: &[Complex<T>],
    lda: isize,
    b: &[Complex<T>],
    ldb: isize,
    beta: Complex<T>,
    c: &mut [Complex<T>],
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
        let mut j = 0;
        while j < n {
            let coorcj = j * ldc;
            let mut i = 0;
            while i < m {
                if beta_is_zero {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else {
                    c[(coorcj + i) as usize] *= beta;
                }
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if side == 'l' || side == 'L' {
        if upper {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let coorcj = j * ldc;
                let mut i = 0;
                while i < m {
                    let coorai = i * lda;
                    let mut tmp = alpha * b[(coorbj + i) as usize];
                    let mut tmp2 = Complex::zero();
                    let mut k = 0;
                    while k < i {
                        c[(coorcj + k) as usize] += tmp * a[(coorai + k) as usize];
                        tmp2 += b[(coorbj + k) as usize] * a[(coorai + k) as usize].conj();
                        k += 1;
                    }
                    tmp = tmp * a[(coorai + i) as usize].re + (alpha * tmp2);
                    if !beta_is_zero {
                        tmp += beta * c[(coorcj + i) as usize];
                    }
                    c[(coorcj + i) as usize] = tmp;
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorbj = j * ldb;
                let coorcj = j * ldc;
                let mut i = m - 1;
                while i >= 0 {
                    let coorai = i * lda;
                    let mut tmp = alpha * b[(coorbj + i) as usize];
                    let mut tmp2 = Complex::zero();
                    let mut k = i + 1;
                    while k < m {
                        c[(coorcj + k) as usize] += tmp * a[(coorai + k) as usize];
                        tmp2 += b[(coorbj + k) as usize] * a[(coorai + k) as usize].conj();
                        k += 1;
                    }
                    tmp = tmp * a[(coorai + i) as usize].re + (alpha * tmp2);
                    if !beta_is_zero {
                        tmp += beta * c[(coorcj + i) as usize];
                    }
                    c[(coorcj + i) as usize] = tmp;
                    i -= 1;
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let cooraj = j * lda;
            let coorbj = j * ldb;
            let coorcj = j * ldc;
            let mut tmp = alpha * a[(cooraj + j) as usize].re;
            let mut i = 0;
            while i < m {
                let mut tmp2 = tmp * b[(coorbj + i) as usize];
                if !beta_is_zero {
                    tmp2 += beta * c[(coorcj + i) as usize];
                }
                c[(coorcj + i) as usize] = tmp2;
                i += 1;
            }
            let mut k = 0;
            while k < j {
                let coorbk = k * ldb;
                let coorak = k * lda;
                if upper {
                    tmp = alpha * a[(cooraj + k) as usize];
                } else {
                    tmp = alpha * a[(coorak + j) as usize].conj();
                };
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] += tmp * b[(coorbk + i) as usize];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1;
            while k < n {
                let coorbk = k * ldb;
                let coorak = k * lda;
                if upper {
                    tmp = alpha * a[(coorak + j) as usize].conj();
                } else {
                    tmp = alpha * a[(cooraj + k) as usize];
                };
                let mut i = 0;
                while i < m {
                    c[(coorcj + i) as usize] += tmp * b[(coorbk + i) as usize];
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
    a: &[Complex<T>],
    lda: isize,
    b: &[Complex<T>],
    ldb: isize,
    beta: T,
    c: &mut [Complex<T>],
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
            let coorcj = j * ldc;
            let f1 = if upper { 0 } else { j };
            let f2 = if upper { j } else { n - 1 };
            if beta.is_zero() {
                let mut i = f1;
                while i <= f2 {
                    c[(coorcj + i) as usize] = Complex::zero();
                    i += 1;
                }
            }
            let start = if upper { 0 } else { j + 1 };
            let stop = if upper { j - 1 } else { n - 1 };
            let mut i = start;
            while i <= stop {
                c[(coorcj + i) as usize] *= beta;
                i += 1;
            }
            let mut ct = c[(coorcj + j) as usize] * beta;
            ct.im = T::zero();
            c[(coorcj + j) as usize] = ct;
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = 0;
            while j < n {
                let coorcj = j * ldc;
                if beta.is_zero() {
                    let mut i = 0;
                    while i <= j {
                        c[(coorcj + i) as usize] = Complex::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < j {
                        c[(coorcj + i) as usize] *= beta;
                        i += 1;
                    }
                    let mut ct = c[(coorcj + j) as usize] * beta;
                    ct.im = T::zero();
                    c[(coorcj + j) as usize] = ct;
                } else {
                    let mut ct = c[(coorcj + j) as usize];
                    ct.im = T::zero();
                    c[(coorcj + j) as usize] = ct;
                }
                let mut l = 0;
                while l < k {
                    let cooral = l * lda;
                    let coorbl = l * ldb;
                    let tmp = b[(coorbl + j) as usize];
                    let tmp2 = a[(cooral + j) as usize];
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        let tmp = alpha * tmp.conj();
                        let tmp2 = (alpha * tmp2).conj();
                        let mut i = 0;
                        while i <= j {
                            c[(coorcj + i) as usize] +=
                                a[(cooral + i) as usize] * tmp + b[(coorbl + i) as usize] * tmp2;
                            i += 1;
                        }
                        let mut ct = c[(coorcj + j) as usize];
                        ct.im = T::zero();
                        c[(coorcj + j) as usize] = ct;
                    }
                    l += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coorcj = j * ldc;
                if beta.is_zero() {
                    let mut i = j;
                    while i < n {
                        c[(coorcj + i) as usize] = Complex::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = j + 1;
                    while i < n {
                        c[(coorcj + i) as usize] *= beta;
                        i += 1;
                    }
                    let mut ct = c[(coorcj + j) as usize] * beta;
                    ct.im = T::zero();
                    c[(coorcj + j) as usize] = ct;
                } else {
                    let mut ct = c[(coorcj + j) as usize];
                    ct.im = T::zero();
                    c[(coorcj + j) as usize] = ct;
                }
                let mut l = 0;
                while l < k {
                    let cooral = l * lda;
                    let coorbl = l * ldb;
                    let mut tmp = b[(coorbl + j) as usize];
                    let mut tmp2 = a[(cooral + j) as usize];
                    if !tmp.is_zero() || !tmp2.is_zero() {
                        tmp = alpha * tmp.conj();
                        tmp2 = (alpha * tmp2).conj();
                        let mut i = j + 1;
                        while i < n {
                            c[(coorcj + i) as usize] +=
                                a[(cooral + i) as usize] * tmp + b[(coorbl + i) as usize] * tmp2;
                            i += 1;
                        }
                        let mut ct = c[(coorcj + j) as usize]
                            + (a[(cooral + j) as usize] * tmp)
                            + (b[(coorbl + j) as usize] * tmp2);
                        ct.im = T::zero();
                        c[(coorcj + j) as usize] = ct
                    }
                    l += 1;
                }
                j += 1;
            }
        }
    } else {
        let mut j = 0;
        while j < n {
            let coorbj = j * ldb;
            let cooraj = j * lda;
            let coorcj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                let coorai = i * lda;
                let coorbi = i * ldb;
                let mut tmp = Complex::zero();
                let mut tmp2 = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coorai + l) as usize].conj() * b[(coorbj + l) as usize];
                    tmp2 += b[(coorbi + l) as usize].conj() * a[(cooraj + l) as usize];
                    l += 1;
                }

                tmp = alpha * tmp + alpha.conj() * tmp2;
                if i == j {
                    tmp.im = T::zero();
                    if beta.is_zero() {
                        c[(coorcj + i) as usize] = tmp;
                    } else {
                        let mut ct = c[(coorcj + i) as usize];
                        ct.im = T::zero();
                        ct *= beta;
                        c[(coorcj + i) as usize] = ct + tmp;
                    }
                } else if beta.is_zero() {
                    c[(coorcj + i) as usize] = tmp;
                } else {
                    c[(coorcj + i) as usize] = c[(coorcj + i) as usize] * beta + tmp;
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
    a: &[Complex<T>],
    lda: isize,
    beta: T,
    c: &mut [Complex<T>],
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
            let coorcj = j * ldc;
            let start = if upper { 0 } else { j };
            let stop = if upper { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                if beta.is_zero() {
                    c[(coorcj + i) as usize] = Complex::zero();
                } else if j == i {
                    // FIXME remove this line if test passes?
                    let mut ct = c[(coorcj + i) as usize] * beta;
                    ct.im = T::zero();
                    c[(coorcj + i) as usize] = ct;
                } else {
                    c[(coorcj + i) as usize] *= beta;
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
            let fe = if upper { j } else { n - 1 };
            let start = if upper { 0 } else { j + 1 };
            let stop = if upper { j - 1 } else { n - 1 };
            let coorcj = j * ldc;
            if beta.is_zero() {
                let mut i = fs;
                while i <= fe {
                    c[(coorcj + i) as usize] = Complex::zero();
                    i += 1;
                }
            } else if !beta.is_one() {
                let mut i = start;
                while i <= stop {
                    c[(coorcj + i) as usize] *= beta;
                    i += 1
                }
                // FIXME only calculate real part.
                let mut ct = c[(coorcj + j) as usize] * beta;
                ct.im = T::zero();
                c[(coorcj + j) as usize] = ct;
            } else {
                let mut ct = c[(coorcj + j) as usize];
                ct.im = T::zero();
                c[(coorcj + j) as usize] = ct;
            }
            let mut l = 0;
            while l < k {
                let cooral = l * lda;
                let mut tmp = a[(cooral + j) as usize];
                if !tmp.is_zero() {
                    tmp = tmp.conj() * alpha;
                    let mut i = start;
                    while i <= stop {
                        c[(coorcj + i) as usize] += tmp * a[(cooral + i) as usize];
                        i += 1;
                    }
                    let mut ct = c[(coorcj + j) as usize] + tmp * a[(cooral + j) as usize];;
                    ct.im = T::zero();
                    c[(coorcj + j) as usize] = ct;
                }
                l += 1;
            }
            j += 1;
        }
    } else if upper {
        let mut j = 0;
        while j < n {
            let cooraj = j * lda;
            let coorcj = j * ldc;
            let mut i = 0;
            while i < j {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coorai + l) as usize].conj() * a[(cooraj + l) as usize];
                    l += 1;
                }
                let mut tmp2 = tmp * alpha;
                if !beta.is_zero() {
                    tmp2 += c[(coorcj + i) as usize] * beta;
                }
                c[(coorcj + i) as usize] = tmp2;
                i += 1;
            }
            let mut rtmp = Complex::zero();
            let mut l = 0;
            while l < k {
                rtmp += a[(cooraj + l) as usize].conj() * a[(cooraj + l) as usize];
                l += 1;
            }
            let mut tmp3: Complex<T> = rtmp * alpha;
            if !beta.is_zero() {
                let re = tmp3.re + beta * c[(coorcj + j) as usize].re;
                tmp3.re = re;
            }
            c[(coorcj + j) as usize] = tmp3;
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let mut rtmp = T::zero();
            let cooraj = j * lda;
            let coorcj = j * ldc;
            let mut l = 0;
            while l < k {
                rtmp += a[(cooraj + l) as usize].norm_sqr();
                l += 1;
            }
            rtmp *= alpha;
            if !beta.is_zero() {
                rtmp += beta * c[(coorcj + j) as usize].re
            }
            c[(coorcj + j) as usize].re = rtmp;
            c[(coorcj + j) as usize].im = T::zero();
            let mut i = j + 1;
            while i < n {
                let coorai = i * lda;
                let mut tmp = Complex::zero();
                let mut l = 0;
                while l < k {
                    tmp += a[(coorai + l) as usize].conj() * a[(cooraj + l) as usize];
                    l += 1;
                }
                tmp *= alpha;
                if !beta.is_zero() {
                    tmp += c[(coorcj + i) as usize] * beta
                }
                c[(coorcj + i) as usize] = tmp;
                i += 1;
            }
            j += 1;
        }
    }
}
