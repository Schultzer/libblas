# libblas
[![documentation](https://docs.rs/libblas/badge.svg)](https://docs.rs/libblas)
[![CircleCI](https://circleci.com/gh/Schultzer/libblas.svg?style=svg)](https://circleci.com/gh/Schultzer/libblas)
[![Build status](https://ci.appveyor.com/api/projects/status/6tywgu4a035iqeqn?svg=true)](https://ci.appveyor.com/project/Schultzer/libblas)


[BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms)for Rust.

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

Conformance testing is done by generating fixtures from [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and run against libblas.
To generate the fixture you need to have [json-fortran](https://github.com/jacobwilliams/json-fortran) installed.

then run `sh ./script/conformance.sh && cargo test --test conformance`

NOTE: only `double precision` and `complex*16` fixtures are genereated.


## Acknowledgement
libblas is implememeted following the [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and [blas-js](https://github.com/R-js/blasjs).


## License
[license]: #license

libblas is primarily distributed under the terms of both the MIT license
and the Apache License (Version 2.0).

See [LICENSE-APACHE](LICENSE-APACHE), [LICENSE-MIT](LICENSE-MIT), and
[COPYRIGHT](COPYRIGHT) for details.