# Student t Log-CDF Gradient

Closed form in the location and the scale, \\-f(q)\\ and \\-z f(q)\\
with \\z = (q-\mu)/\sigma\\; the degrees of freedom are differenced. The
method is
[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
itself, shared with the pseudo-Huber.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of three numeric vectors, `mu`, `sigma` and `nu`, each the
length of `q` recycled against `theta`.

## Details

The derivative of a Student t distribution function with respect to its
degrees of freedom has no elementary form, which is the same obstruction
the skew t meets in its own \\\nu\\ components. One central difference
of the analytic cdf covers it, and only that component pays for the
evaluations.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
degrees of freedom, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.StudentT1Distrib.md)
for the second order;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Examples

``` r
d <- student_t1_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 6)
q <- c(-1, 0.5, 2)

# The location component is exact, the density itself.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# The degrees of freedom are differenced; the component is small and negative
# in the lower tail, heavier tails putting more mass below a low quantile.
distrib_grad_cdf(d, q, th, log = FALSE)$nu
#> [1] -0.0033040362  0.0004428211  0.0039605689
```
