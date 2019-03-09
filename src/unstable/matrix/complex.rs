use num_complex::Complex;
use num_traits::Float;
use std::cmp::{max, min};

pub fn upper_band<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    k: isize,
) -> Vec<Complex<T>> {
    band(mat, rows, cols, 'u', k)
}
pub fn lower_band<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    k: isize,
) -> Vec<Complex<T>> {
    band(mat, rows, cols, 'l', k)
}
pub fn pack_upper<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    k: isize,
) -> Vec<Complex<T>> {
    pack(mat, rows, cols, 'u', k)
}
pub fn pack_lower<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    k: isize,
) -> Vec<Complex<T>> {
    pack(mat, rows, cols, 'l', k)
}

pub fn set_lower<T: Float>(mat: &mut Vec<Complex<T>>, rows: isize, cols: isize, value: T) {
    let mut y = 1;
    while y < cols {
        let coor = (y - 1) * rows - 1;
        let mut i = coor + y + 1;
        while i < coor + rows + 1 {
            mat[i as usize] = Complex {
                re: value,
                im: value,
            };
            i += 1;
        }
        y += 1;
    }
}

pub fn set_upper<T: Float>(mat: &mut Vec<Complex<T>>, rows: isize, cols: isize, value: T) {
    let mut y = 2;
    while y <= cols {
        let coor = (y - 1) * rows - 1;
        let mut i = coor + 1;
        while i < coor + min(y, rows) {
            mat[i as usize] = Complex {
                re: value,
                im: value,
            };
            i += 1;
        }
        y += 1;
    }
}

pub fn slice<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    rowstart: isize,
    rowend: isize,
    colsstart: isize,
    colend: isize,
) -> Vec<Complex<T>> {
    let nrows = (rowend - rowstart) + 1;
    let ncols = (colend - colsstart) + 1;
    let size = (nrows * ncols) as usize;
    let mut data = vec![
        Complex {
            re: T::zero(),
            im: T::zero()
        };
        size
    ];
    let mut j = colsstart;
    while j <= colend {
        let base = (j - colsstart) * nrows - 1;
        let coorj = (j - 1) * rows - 1;
        let mut i = rowstart;
        while i <= rowend {
            let index = (coorj + i) as usize;
            if index == size - 1 {
                // PROPERLY A BUG IN JS. THE FIRST ROUND TURN -1 which is undefined and so it skipz a step.
                // so what happing here is that in the JS code we bind the real with the second imaginary and zero the last one.
                // I have no idea why we do this.
                data[(base + i) as usize] = Complex {
                    re: mat[index].re,
                    im: T::zero(),
                }
            } else {
                data[(base + i) as usize] = Complex {
                    re: mat[index].re,
                    im: mat[index + 1].im,
                }
            }
            i += 1;
        }
        j += 1;
    }
    data
}

fn band<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    uplo: char,
    k: isize,
) -> Vec<Complex<T>> {
    let row_size = k + 1;
    let col_size = cols;
    let size = (row_size * col_size) as usize;
    let mut data = vec![
        Complex {
            re: T::zero(),
            im: T::zero()
        };
        size
    ];
    let mut j = 1;
    while j <= col_size {
        let m = if uplo == 'u' { k + 1 - j } else { 1 - j };
        let coords = (j - 1) * rows - 1;
        let base = (j - 1) * row_size - 1;
        let start = if uplo == 'u' { max(1, j - k) } else { j };
        let stop = if uplo == 'u' { j } else { min(col_size, j + k) };
        let mut i = start;
        while i <= stop {
            data[(base + m + i) as usize] = mat[(coords + i) as usize];
            i += 1;
        }
        j += 1
    }
    data
}

fn pack<T: Float>(
    mat: Vec<Complex<T>>,
    rows: isize,
    cols: isize,
    uplo: char,
    k: isize,
) -> Vec<Complex<T>> {
    let col_size = cols;
    let size = (-(k + 1) * k / 2 + col_size * (k + 1)) as usize;
    let mut data = vec![
        Complex {
            re: T::zero(),
            im: T::zero()
        };
        size
    ];
    let mut cursor = 0;
    let mut j = 1;
    while j <= col_size {
        let coords = (j - 1) * rows - 1;
        let start = if uplo == 'u' { max(1, j - k) } else { j };
        let stop = if uplo == 'u' { j } else { min(col_size, j + k) };
        let mut i = start;
        while i <= stop {
            data[cursor] = mat[(coords + i) as usize];
            i += 1;
            cursor += 1;
        }
        j += 1
    }
    data
}
