# CDF Derivatives of a Sum of Weighted Normal Tails

Evaluates every component of \\\partial^I F\\ of the requested order for
a distribution function of the form \$\$F = c_0 + \sum_k
s_k\\e^{w_k}\\\Phi(x_k).\$\$ Two families in the package have that
shape, the inverse Gaussian and the elastic net, and both reach all four
orders through this one function.

## Usage

``` r
phi_terms_cdf_deriv_k(distrib, q, order, terms)
```

## Arguments

- distrib:

  An object inheriting from `distrib`. Its `params` name and order the
  components.

- q:

  A numeric vector of quantiles.

- order:

  The derivative order, 1 to 4.

- terms:

  A list of terms. Each is a list with `sign` (\\\pm 1\\), `logw` (the
  log weight), `wderiv` (a function of a block of parameter names
  returning that partial of the log weight), `x` (the tail argument) and
  `xderiv` (the same for `x`). A term whose weight is 1 passes a `logw`
  of zero and a `wderiv` returning zero.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## The two sums

The Leibniz rule splits the positions of the multi-index \\I\\ between
the weight and the tail. The weight side is \\\partial^{S} e^{w} = e^{w}
B\_{S}(w)\\, the complete Bell polynomial in the partials of the **log**
weight; the tail side is one Faa di Bruno pass over \\x\\, whose inner
derivatives are the Hermite factors of
[`phi_hermite()`](https://statmodels7.github.io/distributions7/reference/phi_hermite.md)
times \\\varphi(x)\\. Nothing is transcribed: both sums run on the
package's own partition enumeration.

## Why the log weight

The weight is never formed on its own. It is combined with the tail as
`exp(w + pnorm(x, log.p = TRUE))`, because for the inverse Gaussian
\\e^{2/(\phi\mu)}\\ overflows exactly where \\\Phi(b)\\ underflows: at
\\\mu = 0.01\\, \\\phi = 0.1\\ the weight is `Inf` while the product is
about \\3\times10^{-106}\\, and the fourth derivative there comes back
at \\3\times10^{-91}\\.

## Notation

\\F\\ is the distribution function, \\\Phi\\ and \\\varphi\\ the
standard normal distribution and density, \\w\\ a log weight, \\x\\ a
tail argument and \\B_S\\ the complete Bell polynomial.

## See also

[`register_phi_terms_cdf()`](https://statmodels7.github.io/distributions7/reference/register_phi_terms_cdf.md),
which turns a term function into the four methods;
[`phi_hermite()`](https://statmodels7.github.io/distributions7/reference/phi_hermite.md);
[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md)
for the commonest shape of `wderiv` and `xderiv`.
