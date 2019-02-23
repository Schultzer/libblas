# blas-rs
[![documentation](https://docs.rs/blas-rs/badge.svg)](https://docs.rs/blas-rs)
[![CircleCI](https://circleci.com/gh/Schultzer/blas-rs.svg?style=svg)](https://circleci.com/gh/Schultzer/blas-rs)


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


## Acknowledgement
blas-rs is implememeted following the [REFERENCE BLAS Version 3.8.0](http://www.netlib.org/blas/) and [blas-js](https://github.com/R-js/blasjs).


## License
[license]: #license

blas-rs is primarily distributed under the terms of both the MIT license
and the Apache License (Version 2.0).

See [LICENSE-APACHE](LICENSE-APACHE), [LICENSE-MIT](LICENSE-MIT), and
[COPYRIGHT](COPYRIGHT) for details.