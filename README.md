# blas-rs
[![documentation](https://docs.rs/blas-rs/badge.svg)](https://docs.rs/blas-rs)
[![CircleCI](https://circleci.com/gh/Schultzer/blas-rs.svg?style=svg)](https://circleci.com/gh/Schultzer/blas-rs)
[![Build status](https://ci.appveyor.com/api/projects/status/6tywgu4a035iqeqn/branch/master?svg=true)](https://ci.appveyor.com/project/Schultzer/blas-rs/branch/master)


[BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms)for Rust.

## Usage
[usage]: #usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
blas-rs = "0.1"
```

and this to your crate root:

```rust
extern crate blas_rs;
```

## Releases
[releases]: #releases

Release notes are available in [RELEASES.md](RELEASES.md).

## Conformance
[conformance]: #conformance

Conformance testing is done by generating fixtures from [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and run against blas-rs.
To generate the fixture you need to have [json-fortran](https://github.com/jacobwilliams/json-fortran) installed.

then run `sh ./script/conformance.sh && cargo test conformance`

NOTE: only `double precision` and `complex*16` fixtures are genereated.


## Acknowledgement
blas-rs is implememeted following the [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and [blas-js](https://github.com/R-js/blasjs).


## License
[license]: #license

blas-rs is primarily distributed under the terms of both the MIT license
and the Apache License (Version 2.0).

See [LICENSE-APACHE](LICENSE-APACHE), [LICENSE-MIT](LICENSE-MIT), and
[COPYRIGHT](COPYRIGHT) for details.