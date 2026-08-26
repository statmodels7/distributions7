# Derivatives of the Standard Normal Distribution Function

Returns the factor \\h_m\\ in \\\Phi^{(m)}(x) = h_m(x)\\\varphi(x)\\ for
\\m = 1, \ldots, 4\\: \\1\\, \\-x\\, \\x^2 - 1\\ and \\3x - x^3\\. These
are the probabilists' Hermite polynomials up to a sign, and factoring
the density out of them is what keeps the tails evaluable: \\\varphi\\
is supplied by the caller and multiplied in once.

## Usage

``` r
phi_hermite(x, m)
```

## Arguments

- x:

  A numeric vector of arguments.

- m:

  The order, 1 to 4. Any other value returns `NULL`,
  [`switch()`](https://rdrr.io/r/base/switch.html) falling through; no
  caller passes one.

## Value

A numeric vector the length of `x`, the polynomial factor alone.
Multiply by `dnorm(x)` for the derivative itself.

## Notation

\\\Phi\\ is the standard normal distribution function and \\\varphi\\
its density.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md),
the one consumer.

## Examples

``` r
# The first derivative of Phi is the density itself, so h1 is 1.
all.equal(distributions7:::phi_hermite(0.5, 1) * dnorm(0.5), dnorm(0.5))
#> [1] TRUE
```
