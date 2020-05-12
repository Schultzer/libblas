# libblas
[![documentation](https://docs.rs/libblas/badge.svg)](https://docs.rs/libblas)
[![CircleCI](https://circleci.com/gh/Schultzer/libblas.svg?style=svg)](https://circleci.com/gh/Schultzer/libblas)
[![Build status](https://ci.appveyor.com/api/projects/status/6tywgu4a035iqeqn?svg=true)](https://ci.appveyor.com/project/Schultzer/libblas)


[BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) for Rust.

## Why
[why]: #why

 - No special trait representing matrices or vectors
 - Powerful generics
 - Fully tested against [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/)
 - Performance are equal or faster than [Accelerate](https://developer.apple.com/documentation/accelerate) framework
 - No hand tuned assembly

## Examples

```rust
extern crate libblas;
use libblas::level3;

fn main {
  // 3X2 matrix
  let a = vec![1.,2.,3.,4.,5.,6.];
  let b = vec![1.,1.,1.,1.,1.,1.];
  let mut c = vec![2.,1.,3.,1.,4.,1.];
  level3::gemm('n', 'n', 2, 2, 1, 0.3, &a, 3, &b, 3, 1.3, &mut c, 3);
  assert!(c, vec![2.9, 1.9, 3.0, 1.6, 5.8, 1.0])
}
```


## Usage
[usage]: #usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
libblas = "0.1"
```

and this to your crate root:

```rust
extern crate libblas;
```

## Releases
[releases]: #releases

Release notes are available in [RELEASES.md](RELEASES.md).

## Conformance
[conformance]: #conformance

Conformance testing is done by generating fixtures from [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and run against `libblas`.
To generate the fixture you need to have [json-fortran](https://github.com/jacobwilliams/json-fortran) installed.

Run `sh ./script/conformance.sh && cargo test --test conformance`.

NOTE: only `double precision` and `complex*16` fixtures are generated.

## Benchmark
[benchmark]: #benchmark
The benchmark suite does not benchmark against other BLAS implemations. See [libblas-bench](https://github.com/schultzer/libblas-bench).

Run `cargo +nightly bench`

NOTE: remember to have nightly installed `rustup install nightly`.


## Acknowledgement
[acknowledgement]: #acknowledgement

`libblas` is implemented following the [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and [blas-js](https://github.com/R-js/blasjs).


## License
[license]: #license

`libblas` is primarily distributed under the terms of both the MIT license
and the Apache License (Version 2.0).

See [LICENSE-APACHE](LICENSE-APACHE), [LICENSE-MIT](LICENSE-MIT), and
[COPYRIGHT](COPYRIGHT) for details.
