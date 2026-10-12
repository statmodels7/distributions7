# Derivatives of a Discrete Truncated Log-Density Around Its First Point

Computes the parameter derivatives of the log-density of a discrete
parent truncated to \\\\L, L+1, \dots, U\\\\ as \$\$\ell_T(y) =
\\\ell(y) - \ell(L)\\ - \log R, \qquad R = \sum\_{k \ge L} r_k, \quad
r_k = f(k)/f(L).\$\$ The derivatives of \\\ell(y) - \ell(L)\\ and of
\\\log r_k\\ are differences of the parent's own at two points, and
those of \\\log R\\ come from the ratios \\d^B R / R = \sum_k r_k
B_B(\log r_k) / R\\ by the moment-to-cumulant relation of
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md).
Where the retained mass sits at \\L\\ (a zero-truncated count with a
mean near zero), the parent's derivatives and those of \\\log Z\\ are
both of order \\\mu^{-n}\\ and cancel to a number of order one; in this
form the leading power cancels exactly, and the error falls from about
\\\epsilon\mu^{-n}\\ to about \\\epsilon\mu^{1-n}\\.

## Usage

``` r
head_derivs(parent, y, theta, L, U, order, params)
```

## Arguments

- parent:

  A discrete distribution object.

- y:

  A numeric vector of retained observations.

- theta:

  A named list of the parent's parameters, each of length one or of the
  length of `y`.

- L, U:

  The truncation endpoints, `L` finite.

- order:

  The highest derivative order, 1 to 5.

- params:

  The parent's parameter names.

## Value

A function of a block of parameter names returning the derivative of
\\\ell_T\\ with respect to that block, one value per observation.

## Details

The sum over the support starts with eight points and doubles until the
last ratio is below \\e^{-40}\\ in every row, or the support ends, up to
4096 points.

## See also

[`head_override()`](https://statmodels7.github.io/distributions7/reference/head_override.md),
which applies it, and
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md).

## Examples

``` r
p <- poisson_distrib()
h <- distributions7:::head_derivs(p, c(1, 2), list(mu = 1e-6), 1, Inf, 2L, "mu")
h(c("mu", "mu"))
#> [1] -8.333333e-02 -1.000000e+12
```
