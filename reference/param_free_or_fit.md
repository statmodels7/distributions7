# A Free Vector Representing a Matrix, Exactly or Approximately

Returns the free vector of a `parameters7` structure whose matrix is
`m`.
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
answers exactly where the structure can represent the matrix, which for
an unstructured covariance is always. Where it cannot, that call signals
an error and this falls back to least squares on the entries, started
from the zero vector. A compound-symmetric structure asked for an
arbitrary covariance is the ordinary case.

The fallback is a starting value for a starting value and its accuracy
does not matter. It is capped at 200 BFGS iterations and returns the
zero vector if even that fails, so the caller always receives a usable
vector.

## Usage

``` r
param_free_or_fit(s, m)
```

## Arguments

- s:

  A parameters7 structure, supplying `n_free`, `param_free()` and
  `param_value()`.

- m:

  The matrix to represent, of the structure's own dimension.

## Value

An unnamed numeric vector of length `s@n_free`.

## See also

[`mv_moment_start()`](https://statmodels7.github.io/distributions7/reference/mv_moment_start.md),
which supplies `m`;
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
for the exact route.
