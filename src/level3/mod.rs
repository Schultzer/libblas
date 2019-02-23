// use std::arch::x86_64::*;
use num_traits::{Float, NumAssignOps};
use std::cmp::max;

pub mod complex;

/// GEMM  performs one of the matrix-matrix operations
///  C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X   or   op( X ) = X**T,
/// alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
/// This is [SGEMM](http://www.netlib.org/lapack/explore-html/d4/de2/sgemm_8f.html) and [DGEMM](http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html) comined in one function
#[inline]
pub fn gemm<T: Float + NumAssignOps>(
    trans_a: char,
    trans_b: char,
    m: isize,
    n: isize,
    k: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    b: &[T],
    ldb: isize,
    beta: T,
    c: &mut [T],
    ldc: isize,
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
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if k < 0 {
        info = 5;
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

    if m == 0 || n == 0 || (alpha.is_zero() || k == 1) && beta.is_one() {
        return;
    }

    //Refactor #001
    if alpha.is_zero() {
        if beta.is_zero() {
            //Refactor
            let mut j = 0;
            while j < n {
                let coords = j * ldc;
                let mut i = 0;
                while i < m {
                    c[(coords + i) as usize] = T::zero();
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coords = j * ldc;
                let mut i = 0;
                while i < m {
                    c[(coords + i) as usize] *= beta;
                    i += 1;
                }
                j += 1;
            }
        }
        return;
    }

    if not_b {
        if not_a {
            let mut j = 0;
            while j < n {
                let coords_cj = j * ldc;
                let coords_bj = j * ldb;
                if beta.is_zero() {
                    //Refactor
                    let mut i = 0;
                    while i < m {
                        c[(coords_cj + i) as usize] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i < m {
                        c[(coords_cj + i) as usize] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let tmp = alpha * b[(coords_bj + l) as usize];
                    let coords_aj = l * lda;
                    let mut i = 0;
                    while i < m {
                        c[(coords_cj + i) as usize] += tmp * a[(coords_aj + i) as usize];
                        i += 1;
                    }
                    l += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coords_cj = j * ldc;
                let coords_bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let mut tmp = T::zero();
                    let coords = i * lda;
                    let mut l = 0;
                    while l < k {
                        tmp += a[(l + coords) as usize] * b[(coords_bj + l) as usize];
                        l += 1;
                    }
                    if beta.is_zero() {
                        c[(coords_cj + i) as usize] = alpha * tmp;
                    } else {
                        c[(coords_cj + i) as usize] =
                            alpha * tmp + beta * c[(coords_cj + i) as usize];
                    }
                    i += 1;
                }
                j += 1;
            }
        }
    } else if not_a {
        let mut j = 0;
        while j < n {
            let coords_cj = j * ldc;
            if beta.is_zero() {
                //Refactor
                let mut i = 0;
                while i < m {
                    c[(coords_cj + i) as usize] = T::zero();
                    i += 1;
                }
            } else if !beta.is_one() {
                let mut i = 0;
                while i < m {
                    c[(coords_cj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coords_bj = l * ldb;
                let tmp = alpha * b[(coords_bj + j) as usize];
                let coords_aj = l * lda;
                let mut i = 0;
                while i < m {
                    c[(coords_cj + i) as usize] += tmp * a[(coords_aj + i) as usize];
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
                let coords_ai = i * lda;
                let mut l = 0;
                while l < k {
                    let coords_bl = l * ldb;
                    tmp += a[(coords_ai + l) as usize] * b[(coords_bl + j) as usize];
                    l += 1;
                }
                let coords_cj = j * ldc;
                if beta.is_zero() {
                    c[(coords_cj + i) as usize] = alpha * tmp;
                } else {
                    c[(coords_cj + i) as usize] = alpha * tmp + beta * c[(coords_cj + i) as usize];
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
    m: isize,
    n: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    b: &[T],
    ldb: isize,
    beta: T,
    c: &mut [T],
    ldc: isize,
) {
    let nrow_a = if side == 'l' || side == 'L' { m } else { n };
    let mut info = 0;
    if side != 'l' && side != 'L' && side != 'r' && side != 'R' {
        info = 1;
    } else if uplo != 'l' && uplo != 'L' && uplo != 'u' && uplo != 'U' {
        info = 2;
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
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

    //Refactor #001
    if alpha.is_zero() {
        if beta.is_zero() {
            //Refactor
            let mut j = 0;
            while j < n {
                let coords = j * ldc;
                let mut i = 0;
                while i < m {
                    c[(coords + i) as usize] = T::zero();
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coords = j * ldc;
                let mut i = 0;
                while i < m {
                    c[(coords + i) as usize] *= beta;
                    i += 1;
                }
                j += 1;
            }
        }
        return;
    }

    if side == 'l' || side == 'L' {
        if uplo == 'u' || uplo == 'U' {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let coord_cj = j * ldc;
                let mut i = 0;
                while i < m {
                    let coord_ai = i * ldc;
                    let tmp = alpha * b[(coord_bj + i) as usize];
                    let mut tmp2 = T::zero();
                    let mut k = 0;
                    while k < i {
                        c[(coord_cj + k) as usize] += tmp * a[(coord_ai + k) as usize];
                        tmp2 += b[(coord_bj + k) as usize] * a[(coord_ai + k) as usize];
                        k += 1
                    }
                    let mut re = tmp * a[(coord_ai + i) as usize] + alpha * tmp2;
                    if !beta.is_zero() {
                        re += beta * c[(coord_cj + i) as usize];
                    }
                    c[(coord_cj + i) as usize] = re;
                    i += 1
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
                    let coord_ai = i * ldc;
                    let tmp = alpha * b[(coord_bj + i) as usize];
                    let mut tmp2 = T::zero();
                    let mut k = i + 1; // 1
                    while k < m {
                        c[(coord_cj + k) as usize] += tmp * a[(coord_ai + k) as usize];
                        tmp2 += b[(coord_bj + k) as usize] * a[(coord_ai + k) as usize];
                        k += 1
                    }
                    let mut re = tmp * a[(coord_ai + i) as usize] + alpha * tmp2;
                    if !beta.is_zero() {
                        re += beta * c[(coord_cj + i) as usize];
                    }
                    c[(coord_cj + i) as usize] = re;
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
            let tmp = alpha * a[(coord_aj + j) as usize];
            if beta.is_zero() {
                let mut i = 0;
                while i < m {
                    c[(coord_cj + i) as usize] = tmp * b[(coord_bj + i) as usize];
                    i += 1;
                }
            } else {
                let mut i = 0;
                while i < m {
                    c[(coord_cj + i) as usize] =
                        beta * c[(coord_cj + i) as usize] + tmp * b[(coord_bj + i) as usize];
                    i += 1;
                }
            }
            let mut k = 0; //1
            while k < j {
                let coor_ak = k * lda;
                let coor_bk = k * ldb;
                let tmp = if uplo == 'u' || uplo == 'U' {
                    alpha * a[(coord_aj + k) as usize]
                } else {
                    alpha * a[(coor_ak + j) as usize]
                };
                let mut i = 0;
                while i < m {
                    c[(coord_cj + i) as usize] += tmp * b[(coor_bk + i) as usize];
                    i += 1;
                }
                k += 1;
            }
            let mut k = j + 1; //j + 1
            while k < n {
                let coor_ak = k * lda;
                let coor_bk = k * ldb;
                let tmp = if uplo == 'u' || uplo == 'U' {
                    alpha * a[(coor_ak + j) as usize]
                } else {
                    alpha * a[(coord_aj + k) as usize]
                };
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
/// where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
/// and  A and B  are  n by k  matrices  in the  first  case  and  k by n matrices in the second case.
/// This is [SSYR2K](http://www.netlib.org/lapack/explore-html/df/d3d/ssyr2k_8f.html) and [DSYR2K](http://www.netlib.org/lapack/explore-html/d1/dec/dsyr2k_8f.html) comined in one function
#[inline]
pub fn syr2k<T: Float + NumAssignOps>(
    uplo: char,
    trans: char,
    n: isize,
    k: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    b: &[T],
    ldb: isize,
    beta: T,
    c: &mut [T],
    ldc: isize,
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
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
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
            let stop = if uplo == 'u' || uplo == 'U' { j } else { n - 1 };
            let coords = j * ldc;
            if beta.is_zero() {
                while start <= stop {
                    c[(coords + start) as usize] = T::zero();
                    start += 1;
                }
            } else {
                while start <= stop {
                    c[(coords + start) as usize] *= beta;
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
                let coord_cj = j * ldc;
                // Refactor #001
                if beta.is_zero() {
                    // for i in coord_cj + 0c..oord_cj + j + 1 {
                    //   c[i as usize] = T::zero();
                    //   // println!("{:?} {:?}", coord_cj + i, coord_cj + j..coord_cj + j + 1 );
                    // }
                    let mut i = 0;
                    while i <= j {
                        c[(coord_cj + i) as usize] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = 0;
                    while i <= j {
                        c[(coord_cj + i) as usize] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let coord_al = l * lda;
                    let coord_bl = l * ldb;
                    if !a[(coord_al + j) as usize].is_zero()
                        || !b[(coord_bl + j) as usize].is_zero()
                    {
                        let tmp = alpha * b[(coord_bl + j) as usize];
                        let tmp2 = alpha * a[(coord_al + j) as usize];
                        let mut i = 0;
                        while i <= j {
                            c[(coord_cj + i) as usize] += a[(coord_al + i) as usize] * tmp
                                + b[(coord_bl + i) as usize] * tmp2;
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
                let coord_cj = j * ldc;
                // Refactor #001
                if beta.is_zero() {
                    let mut i = j;
                    while i < n {
                        c[(coord_cj + i) as usize] = T::zero();
                        i += 1;
                    }
                } else if !beta.is_one() {
                    let mut i = j;
                    while i < n {
                        c[(coord_cj + i) as usize] *= beta;
                        i += 1;
                    }
                }
                let mut l = 0;
                while l < k {
                    let coord_al = l * lda;
                    let coord_bl = l * ldb;
                    if !a[(coord_al + j) as usize].is_zero()
                        || !b[(coord_bl + j) as usize].is_zero()
                    {
                        let tmp = alpha * b[(coord_bl + j) as usize];
                        let tmp2 = alpha * a[(coord_al + j) as usize];
                        let mut i = j;
                        while i < n {
                            c[(coord_cj + i) as usize] += a[(coord_al + i) as usize] * tmp
                                + b[(coord_bl + i) as usize] * tmp2;
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
            let coord_bj = j * ldb;
            let coord_cj = j * ldc;
            let mut i = 0;
            while i <= j {
                let mut tmp = T::zero();
                let mut tmp2 = T::zero();
                let coord_ai = i * lda;
                let coord_bi = i * ldb;
                let mut l = 0;
                while l < k {
                    tmp += a[(coord_ai + l) as usize] * b[(coord_bj + l) as usize];
                    tmp2 += b[(coord_bi + l) as usize] * b[(coord_bj + l) as usize];
                    l += 1;
                }
                if beta.is_zero() {
                    c[(coord_cj + i) as usize] = alpha * tmp + alpha * tmp2;
                } else {
                    c[(coord_cj + i) as usize] =
                        beta * c[(coord_cj + i) as usize] + alpha * tmp + alpha * tmp2;
                }
                i += 1;
            }
            j += 1;
        }
    } else {
        let mut j = 0;
        while j < n {
            let coord_aj = j * lda;
            let coord_bj = j * ldb;
            let coord_cj = j * ldc;
            let mut i = j;
            while i < n {
                let mut tmp = T::zero();
                let mut tmp2 = T::zero();
                let coord_ai = i * lda;
                let coord_bi = i * ldb;
                let mut l = 0;
                while l < k {
                    tmp += a[(coord_ai + l) as usize] * b[(coord_bj + l) as usize];
                    tmp2 += b[(coord_bi + l) as usize] * a[(coord_aj + l) as usize];
                    l += 1;
                }
                let b = if !beta.is_zero() {
                    beta * c[(coord_cj + i) as usize]
                } else {
                    T::zero()
                };
                c[(coord_cj + i) as usize] = alpha * (tmp + tmp2) + b;
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
    n: isize,
    k: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    beta: T,
    c: &mut [T],
    ldc: isize,
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
    } else if n < 0 {
        info = 3;
    } else if k < 0 {
        info = 4;
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
            let coords = j * ldc;
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j } else { n - 1 };
            let mut i = start;
            while i <= stop {
                c[(coords + i) as usize] = if beta.is_zero() {
                    T::zero()
                } else {
                    beta * c[(coords + i) as usize]
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
            let coord_cj = j * ldc;
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j } else { n - 1 };
            // Refactor #001
            if beta.is_zero() {
                for i in coord_cj + start..=coord_cj + stop {
                    // for i in coord_cj + start..coord_cj + stop + 1 {
                    c[i as usize] = T::zero();
                }
            } else if !beta.is_one() {
                let mut i = start;
                while i <= stop {
                    c[(coord_cj + i) as usize] *= beta;
                    i += 1;
                }
            }
            let mut l = 0;
            while l < k {
                let coord_al = l * lda;
                if !a[(coord_al + j) as usize].is_zero() {
                    let tmp = alpha * a[(coord_al + j) as usize];
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
            let start = if uplo == 'u' || uplo == 'U' { 0 } else { j };
            let stop = if uplo == 'u' || uplo == 'U' { j } else { n - 1 };
            let coord_aj = j * lda;
            let coord_cj = j * ldc;
            let mut i = start;
            while i <= stop {
                let mut tmp = T::zero();
                let coord_ai = i * lda;
                let mut l = 0;
                while l < k {
                    tmp += a[(coord_ai + l) as usize] * a[(coord_aj + l) as usize];
                    l += 1;
                }
                let mut re = alpha * tmp;
                if !beta.is_zero() {
                    re += beta * c[(coord_cj + i) as usize];
                }
                c[(coord_cj + i) as usize] = re;
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
    m: isize,
    n: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    b: &mut [T],
    ldb: isize,
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

    // Refactor #001
    if alpha.is_zero() {
        let mut j = 0;
        while j < n {
            let coords = j * ldb;
            let mut i = 0;
            while i < m {
                b[(coords + i) as usize] = T::zero();
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if lside {
        if trans == 'n' || trans == 'N' {
            if upper {
                let mut j = 0;
                while j < n {
                    let coord_bj = j * ldb;
                    let mut k = 0;
                    while k < m {
                        let coord_ak = k * lda;
                        if !b[(coord_bj + k) as usize].is_zero() {
                            let mut tmp = alpha * b[(coord_bj + k) as usize];
                            let mut i = 0;
                            while i < k {
                                b[(coord_bj + i) as usize] += tmp * a[(coord_ak + i) as usize];
                                i += 1;
                            }
                            if nounit {
                                tmp *= a[(coord_ak + k) as usize]
                            };
                            b[(coord_bj + k) as usize] = tmp;
                        }
                        k += 1;
                    }
                    j += 1;
                }
            } else {
                let mut j = 0;
                while j < n {
                    let coord_bj = j * ldb;
                    let mut k = m - 1;
                    while k >= 0 {
                        let coord_ak = k * lda;
                        if !b[(coord_bj + k) as usize].is_zero() {
                            let tmp = alpha * b[(coord_bj + k) as usize];
                            b[(coord_bj + k) as usize] = tmp;
                            if nounit {
                                b[(coord_bj + k) as usize] *= a[(coord_ak + k) as usize]
                            };
                            let mut i = k + 1;
                            while i < m {
                                b[(coord_bj + i) as usize] += tmp * a[(coord_ak + i) as usize];
                                i += 1;
                            }
                        }
                        k -= 1;
                    }
                    j += 1;
                }
            }
        } else if upper {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let mut i = m - 1;
                while i >= 0 {
                    let mut tmp = b[(coord_bj + i) as usize];
                    let coord_ai = i * lda;
                    if nounit {
                        tmp *= a[(coord_ai + i) as usize]
                    };
                    let mut k = 0;
                    while k < i {
                        tmp += a[(coord_ai + k) as usize] * b[(coord_bj + k) as usize];
                        k += 1;
                    }
                    b[(coord_bj + i) as usize] = alpha * tmp;
                    i -= 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let coord_ai = i * lda;
                    let mut tmp = b[(coord_bj + i) as usize];
                    if nounit {
                        tmp *= a[(coord_ai + i) as usize]
                    }
                    let mut k = i + 1;
                    while k < m {
                        tmp += a[(coord_ai + k) as usize] * b[(coord_bj + k) as usize];
                        k += 1;
                    }
                    b[(coord_bj + i) as usize] = alpha * tmp;
                    i += 1;
                }
                j += 1;
            }
        }
    } else if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = n - 1;
            while j >= 0 {
                let mut tmp = alpha;
                let coord_aj = j * lda;
                let coord_bj = j * ldb;
                if nounit {
                    tmp *= a[(coord_aj + j) as usize]
                };
                let mut i = 0;
                while i < m {
                    b[(coord_bj + i) as usize] *= tmp;
                    i += 1;
                }
                let mut k = 0;
                while k < j {
                    let coord_bk = k * ldb;
                    if !a[(coord_aj + k) as usize].is_zero() {
                        tmp = alpha * a[(coord_aj + k) as usize];
                        let mut i = 0;
                        while i < m {
                            let tmp2 = tmp * b[(coord_bk + i) as usize];
                            b[(coord_bj + i) as usize] += tmp2;
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
                let coord_aj = j * lda;
                let coord_bj = j * ldb;
                let mut tmp = alpha;
                if nounit {
                    tmp *= a[(coord_aj + j) as usize]
                };
                let mut i = 0;
                while i < m {
                    b[(coord_bj + i) as usize] *= tmp;
                    i += 1;
                }
                let mut k = j + 1;
                while k < n {
                    let coord_bk = k * ldb;
                    if !a[(coord_aj + k) as usize].is_zero() {
                        tmp = alpha * a[(coord_aj + k) as usize];
                        let mut i = 0;
                        while i < m {
                            let tmp2 = tmp * b[(coord_bk + i) as usize];
                            b[(coord_bj + i) as usize] += tmp2;
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
            let coord_ak = k * lda;
            let coord_bk = k * ldb;
            let mut j = 0;
            while j < k {
                let coord_bj = j * ldb;
                if !a[(coord_ak + j) as usize].is_zero() {
                    let tmp = alpha * a[(coord_ak + j) as usize];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[(coord_bk + i) as usize];
                        b[(coord_bj + i) as usize] += tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            let mut tmp = alpha;
            if nounit {
                tmp *= a[(coord_ak + k) as usize]
            };
            if tmp != T::one() {
                let mut i = 0;
                while i < m {
                    let tmp2 = tmp * b[(coord_bk + i) as usize];
                    b[(coord_bk + i) as usize] = tmp2;
                    i += 1;
                }
            }
            k += 1;
        }
    } else {
        let mut k = n - 1;
        while k >= 0 {
            let coord_ak = k * lda;
            let coord_bk = k * ldb;
            let mut j = k + 1;
            while j < n {
                let coord_bj = j * ldb;
                if !a[(coord_ak + j) as usize].is_zero() {
                    let tmp = alpha * a[(coord_ak + j) as usize];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[(coord_bk + i) as usize];
                        b[(coord_bj + i) as usize] += tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            let mut tmp = alpha;
            if nounit {
                tmp *= a[(coord_ak + k) as usize]
            };
            if !tmp.is_zero() {
                let mut i = 0;
                while i < m {
                    b[(coord_bk + i) as usize] *= tmp;
                    i += 1
                }
            }
            k -= 1;
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
    m: isize,
    n: isize,
    alpha: T,
    a: &[T],
    lda: isize,
    b: &mut [T],
    ldb: isize,
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

    // Refactor #001
    if alpha.is_zero() {
        let mut j = 0;
        while j < n {
            let coords = j * ldb;
            let mut i = 0;
            while i < m {
                b[(coords + i) as usize] = T::zero();
                i += 1;
            }
            j += 1;
        }
        return;
    }

    if lside {
        if trans == 'n' || trans == 'N' {
            if upper {
                let mut j = 0;
                while j < n {
                    let coord_bj = j * ldb;
                    if !alpha.is_one() {
                        let mut i = 0;
                        while i < m {
                            b[(coord_bj + i) as usize] *= alpha;
                            i += 1;
                        }
                    }
                    let mut k = m - 1;
                    while k >= 0 {
                        let coord_ak = k * lda;
                        if !b[(coord_bj + k) as usize].is_zero() {
                            if nounit {
                                b[(coord_bj + k) as usize] /= a[(coord_ak + k) as usize]
                            };
                            let mut i = 0;
                            while i < k {
                                let tmp = b[(coord_bj + k) as usize] * a[(coord_ak + i) as usize];
                                b[(coord_bj + i) as usize] -= tmp;
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
                    let coord_bj = j * ldb;
                    if !alpha.is_one() {
                        let mut i = 0;
                        while i < m {
                            b[(coord_bj + i) as usize] *= alpha;
                            i += 1;
                        }
                    }
                    let mut k = 0;
                    while k < m {
                        let coord_ak = k * lda;
                        if !b[(coord_bj + k) as usize].is_zero() {
                            if nounit {
                                b[(coord_bj + k) as usize] /= a[(coord_ak + k) as usize]
                            };
                            let mut i = k + 1;
                            while i < m {
                                let tmp = b[(coord_bj + k) as usize] * a[(coord_ak + i) as usize];
                                b[(coord_bj + i) as usize] -= tmp;
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
                let coord_bj = j * ldb;
                let mut i = 0;
                while i < m {
                    let coord_ai = i * lda;
                    let mut tmp = alpha * b[(coord_bj + i) as usize];
                    let mut k = 0;
                    while k < i {
                        tmp -= a[(coord_ai + k) as usize] * b[(coord_bj + k) as usize];
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[(coord_ai + i) as usize]
                    };
                    b[(coord_bj + i) as usize] = tmp;
                    i += 1;
                }
                j += 1;
            }
        } else {
            let mut j = 0;
            while j < n {
                let coord_bj = j * ldb;
                let mut i = m - 1;
                while i >= 0 {
                    let coord_ai = i * lda;
                    let mut tmp = alpha * b[(coord_bj + i) as usize];
                    let mut k = i + 1;
                    while k < m {
                        tmp -= a[(coord_ai + k) as usize] * b[(coord_bj + k) as usize];
                        k += 1;
                    }
                    if nounit {
                        tmp /= a[(coord_ai + i) as usize]
                    };
                    b[(coord_bj + i) as usize] = tmp;
                    i -= 1;
                }
                j += 1;
            }
        }
    } else if trans == 'n' || trans == 'N' {
        if upper {
            let mut j = 0;
            while j < n {
                let coord_aj = j * lda;
                let coord_bj = j * ldb;
                if !alpha.is_one() {
                    let mut i = 0;
                    while i < m {
                        let tmp = alpha * b[(coord_bj + i) as usize];
                        b[(coord_bj + i) as usize] = tmp;
                        i += 1;
                    }
                }
                let mut k = 0;
                while k < j {
                    let coord_bk = k * lda;
                    if !a[(coord_aj + k) as usize].is_zero() {
                        let mut i = 0;
                        while i < m {
                            let tmp = a[(coord_aj + k) as usize] * b[(coord_bk + i) as usize];
                            b[(coord_bj + i) as usize] -= tmp;
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = T::one() / a[(coord_aj + j) as usize];
                    let mut i = 0;
                    while i < m {
                        b[(coord_bj + i) as usize] *= tmp;
                        i += 1;
                    }
                };
                j += 1;
            }
        } else {
            let mut j = n - 1;
            while j >= 0 {
                let coord_aj = j * ldb;
                let coord_bj = j * ldb;
                if !alpha.is_one() {
                    let mut i = 0;
                    while i < m {
                        b[(coord_bj + i) as usize] *= alpha;
                        i += 1;
                    }
                }
                let mut k = j + 1;
                while k < n {
                    let coord_bk = k * ldb;
                    if !a[(coord_aj + k) as usize].is_zero() {
                        let mut i = 0;
                        while i < m {
                            let tmp = a[(coord_aj + k) as usize] * b[(coord_bk + i) as usize];
                            b[(coord_bj + i) as usize] -= tmp;
                            i += 1;
                        }
                    }
                    k += 1;
                }
                if nounit {
                    let tmp = T::one() / a[(coord_aj + j) as usize];
                    let mut i = 0;
                    while i < m {
                        b[(coord_bj + i) as usize] *= tmp;
                        i += 1;
                    }
                };
                j -= 1;
            }
        }
    } else if upper {
        let mut k = n - 1;
        while k >= 0 {
            let coord_ak = k * lda;
            let coord_bk = k * ldb;
            if nounit {
                let tmp = T::one() / a[(coord_ak + k) as usize];
                let mut i = 0;
                while i < m {
                    b[(coord_bk + i) as usize] *= tmp;
                    i += 1;
                }
            };
            let mut j = 0;
            while j < k {
                let coord_bj = j * lda;
                if !a[(coord_ak + j) as usize].is_zero() {
                    let tmp = a[(coord_ak + j) as usize];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[(coord_bk + i) as usize];
                        b[(coord_bj + i) as usize] -= tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha.is_one() {
                let mut i = 0;
                while i < m {
                    b[(coord_bk + i) as usize] *= alpha;
                    i += 1;
                }
            }
            k -= 1;
        }
    } else {
        let mut k = 0;
        while k < n {
            let coord_ak = k * lda;
            let coord_bk = k * ldb;
            if nounit {
                let tmp = T::one() / a[(coord_ak + k) as usize];
                let mut i = 0;
                while i < m {
                    b[(coord_bk + i) as usize] *= tmp;
                    i += 1;
                }
            };
            let mut j = k + 1;
            while j < n {
                let coord_bj = j * ldb;
                if !a[(coord_ak + j) as usize].is_zero() {
                    let tmp = a[(coord_ak + j) as usize];
                    let mut i = 0;
                    while i < m {
                        let tmp2 = tmp * b[(coord_bk + i) as usize];
                        b[(coord_bj + i) as usize] -= tmp2;
                        i += 1;
                    }
                }
                j += 1;
            }
            if !alpha.is_one() {
                let mut i = 0;
                while i < m {
                    b[(coord_bk + i) as usize] *= alpha;
                    i += 1
                }
            }
            k += 1;
        }
    }
}
