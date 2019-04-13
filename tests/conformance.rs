#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;
use libblas::{level1, level2, level3};
use std::fs::File;
use std::io::BufReader;

#[macro_use]
mod utils;
mod case;

#[test]
fn asum() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/asum.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::asum> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::asum(t.n, t.x.as_ptr(), t.incx), t.expect)
    }

    file = File::open("./tests/fixtures/level1/complex/asum.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::asum> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::complex::asum(t.n, t.x.as_ptr(), t.incx), t.expect)
    }
}

#[test]
fn axpy() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/axpy.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::axpy> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level1::axpy(t.n, t.a, t.x.as_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        approximately!(y, t.expect)
    }

    file = File::open("./tests/fixtures/level1/complex/axpy.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::axpy> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level1::complex::axpy(t.n, &t.a, t.x.as_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        capproximately!(y, t.expect)
    }
}

#[test]
fn copy() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/copy.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::copy> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level1::copy(t.n, t.x.as_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        approximately!(y, t.expect)
    }

    file = File::open("./tests/fixtures/level1/complex/copy.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::copy> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level1::complex::copy(t.n, t.x.as_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        capproximately!(y, t.expect)
    }
}

#[test]
fn ddot() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/ddot.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::ddot> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(
            level1::ddot(t.n, t.b, t.x.as_ptr(), t.incx, t.y.as_ptr(), t.incy),
            t.expect
        )
    }
}

#[test]
fn dotc() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/complex/dotc.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::dot> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(
            level1::complex::dotc(t.n, t.x.as_ptr(), t.incx, t.y.as_ptr(), t.incy),
            t.expect
        )
    }
}

#[test]
fn dotu() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/complex/dotu.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::dot> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(
            level1::complex::dotu(t.n, t.x.as_ptr(), t.incx, t.y.as_ptr(), t.incy),
            t.expect
        )
    }
}

#[test]
fn iamax() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/iamax.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::iamax> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::iamax(t.n, t.x.as_ptr(), t.incx), t.expect)
    }

    file = File::open("./tests/fixtures/level1/complex/iamax.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::iamax> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::complex::iamax(t.n, t.x.as_ptr(), t.incx), t.expect)
    }
}

#[test]
fn nrm2() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/nrm2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::nrm2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::nrm2(t.n, t.x.as_ptr(), t.incx), t.expect)
    }

    file = File::open("./tests/fixtures/level1/complex/nrm2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::nrm2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        assert_eq!(level1::complex::nrm2(t.n, t.x.as_ptr(), t.incx), t.expect)
    }
}

#[test]
fn rot() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/rot.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::rot> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        let mut y = t.y;
        level1::rot(
            t.n,
            x.as_mut_ptr(),
            t.incx,
            y.as_mut_ptr(),
            t.incy,
            t.c,
            t.s,
        );
        assert_eq!(x, t.expect.x);
        assert_eq!(y, t.expect.y)
    }

    file = File::open("./tests/fixtures/level1/complex/rot.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::rot> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        let mut y = t.y;
        level1::complex::rot(
            t.n,
            x.as_mut_ptr(),
            t.incx,
            y.as_mut_ptr(),
            t.incy,
            t.c,
            t.s,
        );
        assert_eq!(x, t.expect.x);
        assert_eq!(y, t.expect.y)
    }
}

#[test]
fn rotg() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/rotg.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::rotg> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        let mut b = t.b;
        let mut c = t.c;
        let mut s = t.s;

        level1::rotg(&mut a, &mut b, &mut c, &mut s);
        assert_eq!(a, t.expect.a);
        assert_eq!(b, t.expect.b);
        assert_eq!(c, t.expect.c);
        assert_eq!(s, t.expect.s)
    }

    file = File::open("./tests/fixtures/level1/complex/rotg.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::rotg> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        let mut b = t.b;
        let mut c = t.c;
        let mut s = t.s;
        level1::complex::rotg(&mut a, &mut b, &mut c, &mut s);
        assert_approx_eq_cplx!(a, t.expect.a);
        assert_eq!(b, t.expect.b);
        assert_approx_eq!(c, t.expect.c);
        assert_approx_eq_cplx!(s, t.expect.s)
    }
}

#[test]
fn rotm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/rotm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::rotm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut dx = t.dx;
        let mut dy = t.dy;
        let mut dparam = t.dparam;

        level1::rotm(
            t.n,
            dx.as_mut_ptr(),
            t.incx,
            dy.as_mut_ptr(),
            t.incy,
            dparam.as_mut_ptr(),
        );
        assert_eq!(dx, t.expect.dx);
        assert_eq!(dy, t.expect.dy);
        assert_eq!(dparam, t.expect.dparam);
    }
}

#[test]
fn rotmg() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/rotmg.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::rotmg> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut dd1 = t.dd1;
        let mut dd2 = t.dd2;
        let mut dx1 = t.dx1;
        let mut dy1 = t.dy1;
        let mut dparam = t.dparam;

        level1::rotmg(&mut dd1, &mut dd2, &mut dx1, &mut dy1, dparam.as_mut_ptr());
        assert_approx_eq!(dd1, t.expect.dd1);
        assert_eq!(dd2, t.expect.dd2);
        assert_eq!(dx1, t.expect.dx1);
        assert_eq!(dy1, t.expect.dy1);
        assert_eq!(dparam, t.expect.dparam);
    }
}

#[test]
fn scal() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/scal.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::scal> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level1::scal(t.n, t.a, x.as_mut_ptr(), t.incx);
        assert_eq!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level1/complex/scal.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::scal> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level1::complex::scal(t.n, t.a, x.as_mut_ptr(), t.incx);
        capproximately!(x, t.expect);
    }
}

#[test]
fn sscal() {
    let mut reader;
    let mut file;

    file = File::open("./tests/fixtures/level1/complex/sscal.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::sscal> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level1::complex::sscal(t.n, t.a, x.as_mut_ptr(), t.incx);
        assert_eq!(x, t.expect);
    }
}

#[test]
fn swap() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level1/swap.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::swap> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        let mut y = t.y;
        level1::swap(t.n, x.as_mut_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        assert_eq!(x, t.expect.x);
        assert_eq!(y, t.expect.y);
    }

    file = File::open("./tests/fixtures/level1/complex/swap.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::swap> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        let mut y = t.y;
        level1::complex::swap(t.n, x.as_mut_ptr(), t.incx, y.as_mut_ptr(), t.incy);
        assert_eq!(x, t.expect.x);
        assert_eq!(y, t.expect.y);
    }
}

#[test]
fn gbmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/gbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::gbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::gbmv(
            t.trans, t.m, t.n, t.kl, t.ku, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y,
            t.incy,
        );
        assert_eq!(y, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/gbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::gbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::complex::gbmv(
            t.trans, t.m, t.n, t.kl, t.ku, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y,
            t.incy,
        );
        capproximately!(y, t.expect);
    }
}

#[test]
fn ger() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/ger.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::ger> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        level2::ger(t.m, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut a, t.lda);
        approximately!(a, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/gerc.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::ger> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        println!("{:?}", t.incx);
        level2::complex::gerc(t.m, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut a, t.lda);
        capproximately!(a, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/geru.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::ger> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        level2::complex::geru(t.m, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut a, t.lda);
        capproximately!(a, t.expect);
    }
}

#[test]
fn hbmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/hbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::complex::hbmv(
            t.uplo, t.n, t.k, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        capproximately!(y, t.expect);
    }
}

#[test]
fn hemv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/hemv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hemv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::complex::hemv(
            t.uplo, t.n, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        capproximately!(y, t.expect);
    }
}

#[test]
fn her() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/her.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::her> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        println!("{:?}", t.uplo);
        level2::complex::her(t.uplo, t.n, t.alpha, &t.x, t.incx, &mut a, t.lda);
        capproximately!(a, t.expect);
    }
}

#[test]
fn her2() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/her2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::her2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        level2::complex::her2(
            t.uplo, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut a, t.lda,
        );
        capproximately!(a, t.expect);
    }
}

#[test]
fn hpmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/hpmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hpmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::complex::hpmv(
            t.uplo, t.n, t.alpha, &t.ap, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        capproximately!(y, t.expect);
    }
}

#[test]
fn hpr() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/hpr.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hpr> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut ap = t.ap;
        println!("{:?}", t.uplo);
        level2::complex::hpr(t.uplo, t.n, t.alpha, &t.x, t.incx, &mut ap);
        capproximately!(ap, t.expect);
    }
}

#[test]
fn hpr2() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/complex/hpr2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hpr2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut ap = t.ap;
        level2::complex::hpr2(t.uplo, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut ap);
        capproximately!(ap, t.expect);
    }
}

#[test]
fn sbmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/sbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::sbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::sbmv(
            t.uplo, t.n, t.k, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        approximately!(y, t.expect);
    }
}

#[test]
fn spmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/spmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::spmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::spmv(
            t.uplo, t.n, t.alpha, &t.ap, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        approximately!(y, t.expect);
    }
}

#[test]
fn spr() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/spr.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::spr> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut ap = t.ap;
        level2::spr(t.uplo, t.n, t.alpha, &t.x, t.incx, &mut ap);
        approximately!(ap, t.expect);
    }
}

#[test]
fn spr2() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/spr2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::spr2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut ap = t.ap;
        level2::spr2(t.uplo, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut ap);
        approximately!(ap, t.expect);
    }
}

#[test]
fn symv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/symv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::symv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut y = t.y;
        level2::symv(
            t.uplo, t.n, t.alpha, &t.a, t.lda, &t.x, t.incx, t.beta, &mut y, t.incy,
        );
        approximately!(y, t.expect);
    }
}

#[test]
fn syr() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/syr.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::syr> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        level2::syr(t.uplo, t.n, t.alpha, &t.x, t.incx, &mut a, t.lda);
        approximately!(a, t.expect);
    }
}

#[test]
fn syr2() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/syr2.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::syr2> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut a = t.a;
        level2::syr2(
            t.uplo, t.n, t.alpha, &t.x, t.incx, &t.y, t.incy, &mut a, t.lda,
        );
        approximately!(a, t.expect);
    }
}

#[test]
fn tbmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/tbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::tbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::tbmv(
            t.uplo, t.trans, t.diag, t.n, t.k, &t.a, t.lda, &mut x, t.incx,
        );
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/tbmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::tbmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::tbmv(
            t.uplo, t.trans, t.diag, t.n, t.k, &t.a, t.lda, &mut x, t.incx,
        );
        capproximately!(x, t.expect);
    }
}

#[test]
fn tbsv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/tbsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::tbsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::tbsv(
            t.uplo, t.trans, t.diag, t.n, t.k, &t.a, t.lda, &mut x, t.incx,
        );
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/tbsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::tbsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::tbsv(
            t.uplo, t.trans, t.diag, t.n, t.k, &t.a, t.lda, &mut x, t.incx,
        );
        capproximately!(x, t.expect);
    }
}

#[test]
fn tpmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/tpmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::tpmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::tpmv(t.uplo, t.trans, t.diag, t.n, &t.ap, &mut x, t.incx);
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/tpmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::tpmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::tpmv(t.uplo, t.trans, t.diag, t.n, &t.ap, &mut x, t.incx);
        capproximately!(x, t.expect);
    }
}

#[test]
fn tpsv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/tpsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::tpsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::tpsv(t.uplo, t.trans, t.diag, t.n, &t.ap, &mut x, t.incx);
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/tpsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::tpsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::tpsv(t.uplo, t.trans, t.diag, t.n, &t.ap, &mut x, t.incx);
        capproximately!(x, t.expect);
    }
}

#[test]
fn trmv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/trmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::trmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::trmv(t.uplo, t.trans, t.diag, t.n, &t.a, t.lda, &mut x, t.incx);
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/trmv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::trmv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::trmv(t.uplo, t.trans, t.diag, t.n, &t.a, t.lda, &mut x, t.incx);
        capproximately!(x, t.expect);
    }
}

#[test]
fn trsv() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level2/trsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::trsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::trsv(t.uplo, t.trans, t.diag, t.n, &t.a, t.lda, &mut x, t.incx);
        approximately!(x, t.expect);
    }

    file = File::open("./tests/fixtures/level2/complex/trsv.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::trsv> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut x = t.x;
        level2::complex::trsv(t.uplo, t.trans, t.diag, t.n, &t.a, t.lda, &mut x, t.incx);
        capproximately!(x, t.expect);
    }
}

#[test]
fn gemm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/gemm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::gemm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        println!("{:?}", t.m);
        level3::gemm(
            t.transa, t.transb, t.m, t.n, t.k, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c,
            t.ldc,
        );
        approximately!(c, t.expect, 1E-5);
    }

    file = File::open("./tests/fixtures/level3/complex/gemm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::gemm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::gemm(
            t.transa, t.transb, t.m, t.n, t.k, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c,
            t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn hemm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/complex/hemm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::hemm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::hemm(
            t.side, t.uplo, t.m, t.n, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn herk() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/complex/herk.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::herk> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::herk(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn her2k() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/complex/her2k.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::her2k> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::her2k(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn symm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/symm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::symm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::symm(
            t.side, t.uplo, t.m, t.n, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        approximately!(c, t.expect);
    }

    file = File::open("./tests/fixtures/level3/complex/symm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::symm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::symm(
            t.side, t.uplo, t.m, t.n, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}
#[test]
fn syrk() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/syrk.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::syrk> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::syrk(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, t.beta, &mut c, t.ldc,
        );
        approximately!(c, t.expect);
    }

    file = File::open("./tests/fixtures/level3/complex/syrk.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::syrk> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::syrk(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn syr2k() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/syr2k.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::syr2k> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::syr2k(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        approximately!(c, t.expect);
    }

    file = File::open("./tests/fixtures/level3/complex/syr2k.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::syr2k> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut c = t.c;
        level3::complex::syr2k(
            t.uplo, t.trans, t.n, t.k, t.alpha, &t.a, t.lda, &t.b, t.ldb, t.beta, &mut c, t.ldc,
        );
        capproximately!(c, t.expect);
    }
}

#[test]
fn trmm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/trmm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::trmm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut b = t.b;
        level3::trmm(
            t.side, t.uplo, t.transa, t.diag, t.m, t.n, t.alpha, &t.a, t.lda, &mut b, t.ldb,
        );
        approximately!(b, t.expect);
    }

    file = File::open("./tests/fixtures/level3/complex/trmm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::trmm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut b = t.b;
        level3::complex::trmm(
            t.side, t.uplo, t.transa, t.diag, t.m, t.n, t.alpha, &t.a, t.lda, &mut b, t.ldb,
        );
        capproximately!(b, t.expect);
    }
}

#[test]
fn trsm() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/level3/trsm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::trsm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut b = t.b;
        level3::trsm(
            t.side, t.uplo, t.transa, t.diag, t.m, t.n, t.alpha, &t.a, t.lda, &mut b, t.ldb,
        );
        approximately!(b, t.expect);
    }

    file = File::open("./tests/fixtures/level3/complex/trsm.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::trsm> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        let mut b = t.b;
        level3::complex::trsm(
            t.side, t.uplo, t.transa, t.diag, t.m, t.n, t.alpha, &t.a, t.lda, &mut b, t.ldb,
        );
        capproximately!(b, t.expect);
    }
}

#[test]
fn slice() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/matrix/slice_test.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::matrix> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        capproximately!(t.mat[..6], t.expect);
        capproximately!(t.mat[0..6], t.expect);
    }
}

#[test]
fn pack() {
    let mut reader;
    let mut file;
    file = File::open("./tests/fixtures/matrix/pack_test.json").unwrap();
    reader = BufReader::new(file);
    let tests: Vec<case::complex::matrix> = serde_json::from_reader(reader).unwrap();
    for t in tests {
        //FIXME add more edge cases
        let packed = libblas::unstable::matrix::complex::pack_upper(t.mat.clone(), t.lda, t.n, 0);
        capproximately!(packed, t.expect);
        let packed = libblas::unstable::matrix::complex::pack_lower(t.mat.clone(), t.lda, t.n, 0);
        capproximately!(packed, t.expect);
    }
}
