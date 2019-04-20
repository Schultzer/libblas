pub mod complex;

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
#[serde(default)]
pub struct expect {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub s: f64,
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dd1: f64,
    pub dd2: f64,
    pub dx1: f64,
    pub dy1: f64,
    pub dparam: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct asum {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: f64,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct axpy {
    pub n: isize,
    pub a: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct copy {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct ddot {
    pub n: isize,
    pub b: f32,
    pub x: Vec<f32>,
    pub incx: isize,
    pub y: Vec<f32>,
    pub incy: isize,
    pub expect: f64,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct iamax {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: isize,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct nrm2 {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: f64,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct rot {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub c: f64,
    pub s: f64,
    pub expect: expect,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct rotg {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub s: f64,
    pub expect: expect,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct rotm {
    pub n: isize,
    pub dx: Vec<f64>,
    pub incx: isize,
    pub dy: Vec<f64>,
    pub incy: isize,
    pub dparam: Vec<f64>,
    pub expect: expect,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct rotmg {
    pub dd1: f64,
    pub dd2: f64,
    pub dx1: f64,
    pub dy1: f64,
    pub dparam: Vec<f64>,
    pub expect: expect,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct scal {
    pub n: isize,
    pub a: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct swap {
    pub n: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: expect,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct gbmv {
    pub trans: char,
    pub m: isize,
    pub n: isize,
    pub kl: isize,
    pub ku: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub beta: f64,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct gemv {
    pub trans: char,
    pub m: isize,
    pub n: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub beta: f64,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct ger {
    pub m: isize,
    pub n: isize,
    pub alpha: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct sbmv {
    pub uplo: char,
    pub n: isize,
    pub k: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub beta: f64,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct spmv {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub ap: Vec<f64>,
    pub x: Vec<f64>,
    pub incx: isize,
    pub beta: f64,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct spr {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub ap: Vec<f64>,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct spr2 {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub ap: Vec<f64>,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct symv {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub beta: f64,
    pub y: Vec<f64>,
    pub incy: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct syr {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct syr2 {
    pub uplo: char,
    pub n: isize,
    pub alpha: f64,
    pub x: Vec<f64>,
    pub incx: isize,
    pub y: Vec<f64>,
    pub incy: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct tbmv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub k: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct tbsv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub k: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct tpmv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub ap: Vec<f64>,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct tpsv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub ap: Vec<f64>,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct trmv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}
#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct trsv {
    pub uplo: char,
    pub trans: char,
    pub diag: char,
    pub n: isize,
    pub a: Vec<f64>,
    pub lda: isize,
    pub x: Vec<f64>,
    pub incx: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct gemm {
    pub transa: char,
    pub transb: char,
    pub m: isize,
    pub n: isize,
    pub k: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub b: Vec<f64>,
    pub ldb: isize,
    pub beta: f64,
    pub c: Vec<f64>,
    pub ldc: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct symm {
    pub side: char,
    pub uplo: char,
    pub m: isize,
    pub n: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub b: Vec<f64>,
    pub ldb: isize,
    pub beta: f64,
    pub c: Vec<f64>,
    pub ldc: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct syrk {
    pub uplo: char,
    pub trans: char,
    pub n: isize,
    pub k: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub beta: f64,
    pub c: Vec<f64>,
    pub ldc: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct syr2k {
    pub uplo: char,
    pub trans: char,
    pub n: isize,
    pub k: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub b: Vec<f64>,
    pub ldb: isize,
    pub beta: f64,
    pub c: Vec<f64>,
    pub ldc: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct trmm {
    pub side: char,
    pub uplo: char,
    pub transa: char,
    pub diag: char,
    pub m: isize,
    pub n: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub b: Vec<f64>,
    pub ldb: isize,
    pub expect: Vec<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[allow(non_camel_case_types)]
pub struct trsm {
    pub side: char,
    pub uplo: char,
    pub transa: char,
    pub diag: char,
    pub m: isize,
    pub n: isize,
    pub alpha: f64,
    pub a: Vec<f64>,
    pub lda: isize,
    pub b: Vec<f64>,
    pub ldb: isize,
    pub expect: Vec<f64>,
}
