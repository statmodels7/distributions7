# CDF Derivatives of a Sum of Weighted Normal Tails

Evaluates \\\partial^{I}F\\ for every component of the requested order,
for \\F = c_0 + \sum_k s_k e^{w_k}\Phi(x_k)\\.

## Usage

``` r
phi_terms_cdf_deriv_k(distrib, q, order, terms)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- order:

  The derivative order, 1 to 4.

- terms:

  A list of terms, each a list with `sign`, `logw`, `wderiv`, `x` and
  `xderiv`.

## Value

A named list of derivative components of \\F\\.

## Details

The Leibniz rule splits the positions of \\I\\ between the weight and
the tail, \\\partial^{S}e^{w} = e^{w}B\_{S}(w)\\ being the complete Bell
polynomial in the partials of the log weight and \\\partial^{T}\Phi(x)\\
one Faa di Bruno pass over \\x\\. Nothing is transcribed: both sums run
on the package's own partition enumeration.

## See also

[`bell_f_ratio`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
