# The Data Term of a Multivariate t and Its Derivatives in nu

Computes the \\k\\-th derivative in \\\nu\\ of \\h =
\tfrac{\nu+p}{2}\log(1 + q/\nu)\\, the part of the log-density that
carries the observation, at fixed \\q\\. With \\u = q/\nu\\ and \\a_k(u)
= (1+u)^{-k} - 1\\, for \\k \ge 2\\ \$\$h^{(k)} =
\tfrac{1}{2}(-1)^{k-1}(k-1)!\\\nu^{-(k-1)} \left\\B_k(u) +
\tfrac{p}{\nu}\\a_k(u)\right\\, \qquad B_k(u) = a_k(u) -
\tfrac{k}{k-1}\\a\_{k-1}(u),\$\$ and \\h' = -\tfrac{1}{2}\\D(u) +
(p/\nu)\\u/(1+u)\\\\ with \\D\\ from
[`mvt_D()`](https://statmodels7.github.io/distributions7/reference/mvt_D.md).

## Usage

``` r
mvt_h_d(u, nu, p, k)
```

## Arguments

- u:

  The ratio \\q/\nu\\, a numeric vector.

- nu:

  The degrees of freedom, a single positive number.

- p:

  The dimension.

- k:

  The order of the derivative, 1 to 4.

## Value

A numeric vector as long as `u`.

## Details

The two terms of \\B_k\\ agree to first order in \\u\\, so below \\u =
0.1\\ it is summed as the power series \\\sum\_{m \ge 2}
(-1)^m\\\binom{k+m-1}{m} - \tfrac{k}{k-1} \binom{k+m-2}{m}\\u^m\\, to
\\m = 22\\; above it the direct form loses at most one digit.

## See also

[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
for the part of \\\nu\\ alone; the derivative of the log-density in
\\\nu\\ alone is \\c_p^{(k)} - h^{(k)}\\.

## Examples

``` r
# Against the direct form at nu = 1e9, where its terms cancel.
nu <- 1e9; q <- 1.25; p <- 2
h <- function(v) (v + p) / 2 * log1p(q / v)
distributions7:::mvt_h_d(q / nu, nu, p, 2L)
#> [1] 1.71875e-27
(h(nu * (1 + 1e-4)) - 2 * h(nu) + h(nu * (1 - 1e-4))) / (nu * 1e-4)^2
#> [1] -1.110223e-26
```
