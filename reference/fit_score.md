# Summed Score on the Link Scale

Returns \\\partial \ell / \partial \eta\\ for the whole sample: the
package's per-observation gradient on the link scale, summed over
observations, as a plain numeric vector in parameter order. At the
maximum likelihood estimate every entry is zero to rounding, and a
Gaussian fitted to 200 draws answers about \\10^{-15}\\ there against
204 and 752 one unit away.

## Usage

``` r
fit_score(distrib, y, theta, threads = 1L)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations, or the response matrix of a
  multivariate family.

- theta:

  A named list of parameters on the parameter scale, aligned to
  `distrib@params`.

- threads:

  How many threads the family's compiled kernels may use, as an integer
  count. Defaults to `1L`, which takes the sequential path. The sum does
  not depend on the count.

## Value

An unnamed numeric vector of length `length(distrib@params)`, in the
order `distrib@params` gives.

## See also

[`fit_hess_matrix()`](https://statmodels7.github.io/distributions7/reference/fit_hess_matrix.md)
for the second-order counterpart;
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
which supplies the per-observation components;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
the caller.
