# Gaussian Log-CDF Gradient

Computes the gradient in closed form, from the location-scale structure:
\\\partial F/\partial\mu = -f(q)\\ and \\\partial F/\partial\sigma = -z
f(q)\\ with \\z = (q-\mu)/\sigma\\. The method is
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
itself, shared with the logistic, the Laplace and the Gumbel, so the
method uses nothing particular to the Gaussian beyond its density.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma` (positive); each is a
  numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail. With `log = FALSE`, `FALSE` flips the sign of
  every component; with `log = TRUE` it returns the derivatives of
  \\\log S\\, where \\S = 1 - F\\, which are not the negatives of the
  lower-tail ones.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability. The ratio
  \\f/P\\ is formed as \\\exp(\log f - \log P)\\, so the result stays
  finite far into either tail, where \\f\\ and \\P\\ underflow.

## Value

A named list of two numeric vectors, `mu` and `sigma`, in which each
vector has the length of `q` recycled against `theta`.

## Details

A censored Gaussian likelihood is cheap for this reason: a
right-censored observation contributes \\\log S(q)\\, whose score in
\\\mu\\ is \\f(q)/S(q)\\, the inverse Mills ratio of \\z\\ divided by
\\\sigma\\, and the whole of it comes from the density that the family
already computes.

## Notation

\\\mu\\ is the mean, \\\sigma \> 0\\ the standard deviation, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\F\\ the distribution function.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Gaussian1Distrib.md)
for the second order;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# On the natural scale the mean component is minus the density.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# On the log scale it is -f/F, and the upper tail flips to +f/S.
distrib_grad_cdf(d, q, th)$mu
#> [1] -1.3268963 -0.5790812 -0.1322307
distrib_grad_cdf(d, q, th, lower.tail = FALSE)$mu
#> [1] 0.2148057 0.7557727 1.5567503
```
