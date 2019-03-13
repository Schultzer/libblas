// use std::arch::x86_64::*;
use num_traits::{Float, NumAssignOps};
use std::cmp::max;

pub mod complex;

fn multiply<T: Float + NumAssignOps>(left: &mut [T], right: T, ld: usize, n: usize, m: usize) {
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

fn zero<T: Float + NumAssignOps>(left: &mut [T], ld: usize, n: usize, m: usize) {
    let mut j = 0;
    while j < n {
        let coords = j * ld;
        let mut i = 0;
        while i < m {
            left[coords + i] = T::zero();
            i += 1;
        }
        j += 1;
    }
}

/// GEMM  performs one of the matrix-matrix operations
///  C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X   or   op( X ) = X**T,
/// alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
/// This is [SGEMM](http://www.netlib.org/lapack/explore-html/d4/de2/sgemm_8f.html) and [DGEMM](http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html) comined in one function
#[inline]
pub fn gemm<T: Float + NumAssignOps>(
    trans_a: char,
    trans_b: char,
    m: usize,
    n: usize,
    k: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    b: &[T],
    ldb: usize,
    beta: T,
    c: &mut [T],
    ldc: usize,
) {
    let not_a = trans_a == 'n' || trans_a == 'N';
    let not_b = trans_b == 'n' || trans_b == 'N';
    let nrow_a = if not_a { m } else { k };
    let nrow_b = if not_a { k } else { n };

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
    } else if lda < max(1, nrow_a) {
        info = 8;
    } else if ldb < max(1, nrow_b) {
        info = 10;
    } else if ldc < max(1, m) {
        info = 13;
    }
    if info != 0 {
        panic!("gemm {}", info);
    }

    if m == 0 || n == 0 || (alpha.is_zero() || k == 0) && beta.is_one() {
        return;
    }

    if alpha.is_zero() {
        multiply(c, beta, ldc, n, m);
        return;
    }

    if not_b {
        if not_a {
            let mut j = 0;
            while j < n {
                let cj = j * ldc;
                let bj = j * ldb;
                if beta.is_zero() {
                    //Refactor
                    let mut i = 0;
                    while i < m {
                        c[cj + i] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < m {
                        c[cj + i] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let tmp = alpha * b[bj + l];
                    let aj = l * lda;
                    let mut i = 0;
                    while i < m {
                        c[cj + i] += tmp * a[aj + i];
                        i += 1;
                    }
                    l += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let cj = j * ldc;
                let bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let mut tmp = T::zero();
                    let coords = i * lda;
                    let mut l = 0;
                    while l < k {
                        tmp += a[l + coords] * b[bj + l];
                        l += 1;
                    }
                    if beta.is_zero() {
                        c[cj + i] = alpha * tmp;
                    } else {
                        c[cj + i] = alpha * tmp + beta * c[cj + i];
                    }
                    i += 1;
                }
                j += 1;
            }
        }
    } else if not_a {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            if beta.is_zero() {
                //Refactor
                let mut i = 0;
                while i < m {
                    c[cj + i] = T::zero();
                    i += 1;
                }
            } else if !beta.is_one() {
                let mut i = 0;
                while i < m {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let bj = l * ldb;
                let tmp = alpha * b[bj + j];
                let aj = l * lda;
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * a[aj + i];
                    i += 1;
                }
                l += 1;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let mut i = 0;
            while i < m {
                let mut tmp = T::zero();
                let ai = i * lda;
                let mut l = 0;
                while l < k {
                    let bl = l * ldb;
                    tmp += a[ai + l] * b[bl + j];
                    l += 1;
                }
                let cj = j * ldc;
                if beta.is_zero() {
                    c[cj + i] = alpha * tmp;
                } else {
                    c[cj + i] = alpha * tmp + beta * c[cj + i];
                }
                i += 1;
            }
            j += 1;
        }
    }
}

/// SYMM  performs one of the matrix-matrix operations
/// C := alpha*A*B + beta*C, or  C := alpha*B*A + beta*C,
/// where alpha and beta are scalars,  A is a symmetric matrix and  B and C are  m by n matrices.
/// This is [SSYMM](http://www.netlib.org/lapack/explore-html/d7/d42/ssymm_8f.html) and [DSYMM](http://www.netlib.org/lapack/explore-html/d8/db0/dsymm_8f.html) comined in one function
#[inline]
pub fn symm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    m: usize,
    n: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    b: &[T],
    ldb: usize,
    beta: T,
    c: &mut [T],
    ldc: usize,
) {
    let nrow_a = if side == 'l' || side == 'L' { m } else { n };
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if lda < max(1, nrow_a) {
        info = 7;
    } else if ldb < max(1, m) {
        info = 9;
    } else if ldc < max(1, m) {
        info = 12;
    }
    if info != 0 {
        panic!("symm {}", info);
    }

    if m == 0 || n == 0 || (alpha.is_zero() && beta.is_one()) {
        return;
    }

    if alpha.is_zero() {
        multiply(c, beta, ldc, n, m);
        return;
    }

    if side == 'l' || side == 'L' {
        if uplo == 'u' || uplo == 'U' {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let cj = j * ldc;
                let mut i = 0;
                while i < m {
                    let ai = i * ldc;
                    let tmp = alpha * b[bj + i];
                    let mut tmp2 = T::zero();
                    let mut k = 0;
                    while k < i {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k];
                        k += 1
                    }
                    let mut re = tmp * a[ai + i] + alpha * tmp2;
                    if !beta.is_zero() {
                        re += beta * c[cj + i];
                    }
                    c[cj + i] = re;
                    i += 1
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
                    let ai = i * ldc;
                    let tmp = alpha * b[bj + i];
                    let mut tmp2 = T::zero();
                    let mut k = i + 1;
                    while k < m {
                        c[cj + k] += tmp * a[ai + k];
                        tmp2 += b[bj + k] * a[ai + k];
                        k += 1
                    }
                    let mut re = tmp * a[ai + i] + alpha * tmp2;
                    if !beta.is_zero() {
                        re += beta * c[cj + i];
                    }
                    c[cj + i] = re;
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
            let tmp = alpha * a[aj + j];
            if beta.is_zero() {
                let mut i = 0;
                while i < m {
                    c[cj + i] = tmp * b[bj + i];
                    i += 1;
                }
            } else {
                let mut i = 0;
                while i < m {
                    c[cj + i] = beta * c[cj + i] + tmp * b[bj + i];
                    i += 1;
                }
            }
            let mut k = 0; //1
            while k < j {
                let ak = k * lda;
                let bk = k * ldb;
                let tmp = if uplo == 'u' || uplo == 'U' {
                    alpha * a[aj + k]
                } else {
                    alpha * a[ak + j]
                };
                let mut i = 0;
                while i < m {
                    c[cj + i] += tmp * b[bk + i];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1; //j + 1
            while k < n {
                let ak = k * lda;
                let bk = k * ldb;
                let tmp = if uplo == 'u' || uplo == 'U' {
                    alpha * a[ak + j]
                } else {
                    alpha * a[aj + k]
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

/// SYR2K  performs one of the symmetric rank 2k operations
/// C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C,
/// where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
/// and  A and B  are  n by k  matrices  in the  first  case  and  k by n matrices in the second case.
/// This is [SSYR2K](http://www.netlib.org/lapack/explore-html/df/d3d/ssyr2k_8f.html) and [DSYR2K](http://www.netlib.org/lapack/explore-html/d1/dec/dsyr2k_8f.html) comined in one function
#[inline]
pub fn syr2k<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: usize,
    k: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    b: &[T],
    ldb: usize,
    beta: T,
    c: &mut [T],
    ldc: usize,
) {
    let nrow_a = if trans == 'n' || trans == 'N' { n } else { k };
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
    } else if lda < max(1, nrow_a) {
        info = 7;
    } else if ldb < max(1, nrow_a) {
        info = 9;
    } else if ldc < max(1, n) {
        info = 12;
    }

    if info != 0 {
        panic!("syr2k {}", info);
    }

    if n == 0 || ((alpha.is_zero() || k == 0) && beta.is_one()) {
        return;
    }

    //Refactor 001
    if alpha.is_zero() {
        let mut j = 0;
        while j < n {
            let mut start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j + 1 } else { n };
            let cj = j * ldc;
            if beta.is_zero() {
                while start < stop {
                    c[cj + start] = T::zero();
                    start += 1;
                }
            } else {
                while start < stop {
                    c[cj + start] *= beta;
                    start += 1;
                }
            }
            j += 1;
        }
        return;
    }

    if trans == 'n' || trans == 'N' {
        if uplo == 'u' || uplo == 'U' {
            let mut j = 0;
            while j < n {
                let cj = j * ldc;
                // Refactor #001
                if beta.is_zero() {
                    let mut i = 0;
                    while i < j + 1 {
                        c[cj + i] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < j + 1 {
                        c[cj + i] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    if !a[al + j].is_zero() || !b[bl + j].is_zero() {
                        let tmp = alpha * b[bl + j];
                        let tmp2 = alpha * a[al + j];
                        let mut i = 0;
                        while i < j + 1 {
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
                let cj = j * ldc;
                // Refactor #001
                if beta.is_zero() {
                    let mut i = j;
                    while i < n {
                        c[cj + i] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = j;
                    while i < n {
                        c[cj + i] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let al = l * lda;
                    let bl = l * ldb;
                    if !a[al + j].is_zero() || !b[bl + j].is_zero() {
                        let tmp = alpha * b[bl + j];
                        let tmp2 = alpha * a[al + j];
                        let mut i = j;
                        while i < n {
                            c[cj + i] += a[al + i] * tmp + b[bl + i] * tmp2;
                            i += 1;
                        }
                    }
                    l += 1;
                }
                j += 1;
            }
        }
    } else if uplo == 'u' || uplo == 'U' {
        let mut j = 0;
        while j < n {
            let bj = j * ldb;
            let cj = j * ldc;
            let mut i = 0;
            while i < j + 1 {
                let mut tmp = T::zero();
                let mut tmp2 = T::zero();
                let ai = i * lda;
                let bi = i * ldb;
                let mut l = 0;
                while l < k {
                    tmp += a[ai + l] * b[bj + l];
                    tmp2 += b[bi + l] * b[bj + l];
                    l += 1;
                }
                if beta.is_zero() {
                    c[cj + i] = alpha * tmp + alpha * tmp2;
                } else {
                    c[cj + i] = beta * c[cj + i] + alpha * tmp + alpha * tmp2;
                }
                i += 1;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let aj = j * lda;
            let bj = j * ldb;
            let cj = j * ldc;
            let mut i = j;
            while i < n {
                let mut tmp = T::zero();
                let mut tmp2 = T::zero();
                let ai = i * lda;
                let bi = i * ldb;
                let mut l = 0;
                while l < k {
                    tmp += a[ai + l] * b[bj + l];
                    tmp2 += b[bi + l] * a[aj + l];
                    l += 1;
                }
                let b = if !beta.is_zero() {
                    beta * c[cj + i]
                } else {
                    T::zero()
                };
                c[cj + i] = alpha * (tmp + tmp2) + b;
                i += 1;
            }
            j += 1;
        }
    }
}

/// SYRK  performs one of the symmetric rank k operations
/// C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C,
/// where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
/// and  A  is an  n by k  matrix in the first case and a  k by n  matrix in the second case.
/// This is [SSYRK](http://www.netlib.org/lapack/explore-html/d0/d40/ssyrk_8f.html) and [DSYRK](http://www.netlib.org/lapack/explore-html/dc/d05/dsyrk_8f.html) comined in one function
#[inline]
pub fn syrk<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: usize,
    k: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    beta: T,
    c: &mut [T],
    ldc: usize,
) {
    let nrow_a = if trans == 'n' || trans == 'N' { n } else { k };
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
    } else if lda < max(1, nrow_a) {
        info = 7;
    } else if ldc < max(1, n) {
        info = 10;
    }

    if info != 0 {
        panic!("syrk {}", info);
    }

    if n == 0 || ((alpha.is_zero() || k == 0) && beta.is_one()) {
        return;
    }

    //Refactor 001
    if alpha.is_zero() {
        let mut j = 0;
        while j < n {
            let cj = j * ldc;
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j + 1 } else { n };
            let mut i = start;
            while i < stop {
                c[cj + i] = if beta.is_zero() {
                    T::zero()
                } else {
                    beta * c[cj + i]
                };
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
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j + 1 } else { n };
            // Refactor #001
            if beta.is_zero() {
                for i in cj + start..cj + stop {
                    // for i in cj + start..cj + stop + 1 {
                    c[i] = T::zero();
                }
            } else if !beta.is_one() {
                let mut i = start;
                while i < stop {
                    c[cj + i] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let al = l * lda;
                if !a[al + j].is_zero() {
                    let tmp = alpha * a[al + j];
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
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j + 1 } else { n };
            let aj = j * lda;
            let cj = j * ldc;
            let mut i = start;
            while i < stop {
                let mut tmp = T::zero();
                let ai = i * lda;
                let mut l = 0;
                while l < k {
                    tmp += a[ai + l] * a[aj + l];
                    l += 1;
                }
                let mut re = alpha * tmp;
                if !beta.is_zero() {
                    re += beta * c[cj + i];
                }
                c[cj + i] = re;
                i += 1;
            }
            j += 1;
        }
    }
}

/// TRMM  performs one of the matrix-matrix operations
/// B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
/// where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
///  op( A ) = A   or   op( A ) = A**T.
/// This is [STRMM](http://www.netlib.org/lapack/explore-html/df/d01/strmm_8f.html) and [DTRMM](http://www.netlib.org/lapack/explore-html/dd/d19/dtrmm_8f.html) comined in one function
#[inline]
pub fn trmm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    trans: char,
    diag: char,
    m: usize,
    n: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    b: &mut [T],
    ldb: usize,
) {
    let lside = side == 'l' || side == 'L';
    let nrowa = if lside { m } else { n };
    let nounit = diag == 'n' || diag == 'N';
    let upper = uplo == 'u' || uplo == 'U';
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

    if alpha.is_zero() {
        zero(b, ldb, n, m);
        return;
    }

    if lside {
        if trans == 'n' || trans == 'N' {
            if upper {
                let mut j = 0;
                while j < n {
                    let bj = j * ldb;
                    let mut k = 0;
                    while k < m {
                        let ak = k * lda;
                        if !b[bj + k].is_zero() {
                            let mut tmp = alpha * b[bj + k];
                            let mut i = 0;
                            while i < k {
                                b[bj + i] += tmp * a[ak + i];
                                i += 1;
                            }
                            if nounit {
                                tmp *= a[ak + k]
                            };
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
                        let ak = k * lda;
                        if !b[bj + k].is_zero() {
                            let tmp = alpha * b[bj + k];
                            b[bj + k] = tmp;
                            if nounit {
                                b[bj + k] *= a[ak + k]
                            };
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
        } else if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = m;
                while i >= 1 {
                    i -= 1;
                    let mut tmp = b[bj + i];
                    let ai = i * lda;
                    if nounit {
                        tmp *= a[ai + i]
                    };
                    let mut k = 0;
                    while k < i {
                        tmp += a[ai + k] * b[bj + k];
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
                        tmp *= a[ai + i]
                    }
                    let mut k = i + 1;
                    while k < m {
                        tmp += a[ai + k] * b[bj + k];
                        k += 1;
                    }
                    b[bj + i] = alpha * tmp;
                    i += 1;
                }
                j += 1;
            }
        }
    } else if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let mut tmp = alpha;
                let aj = j * lda;
                let bj = j * ldb;
                if nounit {
                    tmp *= a[aj + j]
                };
                let mut i = 0;
                while i < m {
                    b[bj + i] *= tmp;
                    i += 1;
                }
                let mut k = 0;
                while k < j {
                    let bk = k * ldb;
                    if !a[aj + k].is_zero() {
                        tmp = alpha * a[aj + k];
                        let mut i = 0;
                        while i < m {
                            let tmp2 = tmp * b[bk + i];
                            b[bj + i] += tmp2;
                            i += 1;
                        }
                    }
                    k += 1;
                }
            }
        } else {
            let mut j = 0;
            while j < n {
                let aj = j * lda;
                let bj = j * ldb;
                let mut tmp = alpha;
                if nounit {
                    tmp *= a[aj + j]
                };
                let mut i = 0;
                while i < m {
                    b[bj + i] *= tmp;
                    i += 1;
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    if !a[aj + k].is_zero() {
                        tmp = alpha * a[aj + k];
                        let mut i = 0;
                        while i < m {
                            let tmp2 = tmp * b[bk + i];
                            b[bj + i] += tmp2;
                            i += 1;
                        }
                    }
                    k += 1;
                }
                j += 1;
            }
        }
    } else if upper {
        let mut k = 0;
        while k < n {
            let ak = k * lda;
            let bk = k * ldb;
            let mut j = 0;
            while j < k {
                let bj = j * ldb;
                if !a[ak + j].is_zero() {
                    let tmp = alpha * a[ak + j];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[bk + i];
                        b[bj + i] += tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            let mut tmp = alpha;
            if nounit {
                tmp *= a[ak + k]
            };
            if tmp != T::one() {
                let mut i = 0;
                while i < m {
                    let tmp2 = tmp * b[bk + i];
                    b[bk + i] = tmp2;
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
                let bj = j * ldb;
                if !a[ak + j].is_zero() {
                    let tmp = alpha * a[ak + j];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[bk + i];
                        b[bj + i] += tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            let mut tmp = alpha;
            if nounit {
                tmp *= a[ak + k]
            };
            if !tmp.is_zero() {
                let mut i = 0;
                while i < m {
                    b[bk + i] *= tmp;
                    i += 1
                }
            }
        }
    }
}

/// TRSM  solves one of the matrix equations
/// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
/// where alpha is a scalar, X and B are m by n matrices, A is a unit, or non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
/// op( A ) = A   or   op( A ) = A**T. The matrix X is overwritten on B.
/// This is [STRSM](http://www.netlib.org/lapack/explore-html/d2/d8b/strsm_8f.html) and [DTRSM](http://www.netlib.org/lapack/explore-html/de/da7/dtrsm_8f.html) comined in one function
#[inline]
pub fn trsm<T: Float + NumAssignOps>(
    side: char,
    uplo: char,
    trans: char,
    diag: char,
    m: usize,
    n: usize,
    alpha: T,
    a: &[T],
    lda: usize,
    b: &mut [T],
    ldb: usize,
) {
    let lside = side == 'l' || side == 'L';
    let nrowa = if lside { m } else { n };
    let nounit = diag == 'n' || diag == 'N';
    let upper = uplo == 'u' || uplo == 'U';
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

    if alpha.is_zero() {
        zero(b, ldb, n, m);
        return;
    }

    if lside {
        if trans == 'n' || trans == 'N' {
            if upper {
                let mut j = 0;
                while j < n {
                    let bj = j * ldb;
                    if !alpha.is_one() {
                        let mut i = 0;
                        while i < m {
                            b[bj + i] *= alpha;
                            i += 1;
                        }
                    }
                    let mut k = m;
                    while k >= 1 {
                        k -= 1;
                        let ak = k * lda;
                        if !b[bj + k].is_zero() {
                            if nounit {
                                b[bj + k] /= a[ak + k]
                            };
                            let mut i = 0;
                            while i < k {
                                let tmp = b[bj + k] * a[ak + i];
                                b[bj + i] -= tmp;
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
                    if !alpha.is_one() {
                        let mut i = 0;
                        while i < m {
                            b[bj + i] *= alpha;
                            i += 1;
                        }
                    }
                    let mut k = 0;
                    while k < m {
                        let ak = k * lda;
                        if !b[bj + k].is_zero() {
                            if nounit {
                                b[bj + k] /= a[ak + k]
                            };
                            let mut i = k + 1;
                            while i < m {
                                let tmp = b[bj + k] * a[ak + i];
                                b[bj + i] -= tmp;
                                i += 1;
                            }
                        }
                        k += 1;
                    }
                    j += 1;
                }
            }
        } else if upper {
            let mut j = 0;
            while j < n {
                let bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let ai = i * lda;
                    let mut tmp = alpha * b[bj + i];
                    let mut k = 0;
                    while k < i {
                        tmp -= a[ai + k] * b[bj + k];
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[ai + i]
                    };
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
                        tmp -= a[ai + k] * b[bj + k];
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[ai + i]
                    };
                    b[bj + i] = tmp;
                }
                j += 1;
            }
        }
    } else if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = 0;
            while j < n {
                let aj = j * lda;
                let bj = j * ldb;
                if !alpha.is_one() {
                    let mut i = 0;
                    while i < m {
                        let tmp = alpha * b[bj + i];
                        b[bj + i] = tmp;
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < j {
                    let bk = k * lda;
                    if !a[aj + k].is_zero() {
                        let mut i = 0;
                        while i < m {
                            let tmp = a[aj + k] * b[bk + i];
                            b[bj + i] -= tmp;
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = T::one() / a[aj + j];
                    let mut i = 0;
                    while i < m {
                        b[bj + i] *= tmp;
                        i += 1;
                    }
                };
                j += 1;
            }
        } else {
            let mut j = n;
            while j >= 1 {
                j -= 1;
                let aj = j * ldb;
                let bj = j * ldb;
                if !alpha.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[bj + i] *= alpha;
                        i += 1;
                    }
                }
                let mut k = j + 1;
                while k < n {
                    let bk = k * ldb;
                    if !a[aj + k].is_zero() {
                        let mut i = 0;
                        while i < m {
                            let tmp = a[aj + k] * b[bk + i];
                            b[bj + i] -= tmp;
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = T::one() / a[aj + j];
                    let mut i = 0;
                    while i < m {
                        b[bj + i] *= tmp;
                        i += 1;
                    }
                };
            }
        }
    } else if upper {
        let mut k = n;
        while k >= 1 {
            k -= 1;
            let ak = k * lda;
            let bk = k * ldb;
            if nounit {
                let tmp = T::one() / a[ak + k];
                let mut i = 0;
                while i < m {
                    b[bk + i] *= tmp;
                    i += 1;
                }
            };
            let mut j = 0;
            while j < k {
                let bj = j * lda;
                if !a[ak + j].is_zero() {
                    let tmp = a[ak + j];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[bk + i];
                        b[bj + i] -= tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha.is_one() {
                let mut i = 0;
                while i < m {
                    b[bk + i] *= alpha;
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
                let tmp = T::one() / a[ak + k];
                let mut i = 0;
                while i < m {
                    b[bk + i] *= tmp;
                    i += 1;
                }
            };
            let mut j = k + 1;
            while j < n {
                let bj = j * ldb;
                if !a[ak + j].is_zero() {
                    let tmp = a[ak + j];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[bk + i];
                        b[bj + i] -= tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha.is_one() {
                let mut i = 0;
                while i < m {
                    b[bk + i] *= alpha;
                    i += 1
                }
            }
            k += 1;
        }
    }
}
