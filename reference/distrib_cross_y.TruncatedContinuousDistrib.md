# Mixed Derivatives of a Truncated Distribution

The parent's mixed derivatives, unchanged. Truncation replaces the
parent's log-density by \\\ell(y;\theta) - \log Z(\theta)\\, and the
normalizing constant does not depend on \\y\\, so it vanishes from any
derivative taking one derivative in the response. The equality is exact,
not approximate, and holds at every observation inside the support.

This is the one derivative of a truncated family that costs nothing. The
gradient and the Hessian both need \\\log Z\\ and its derivatives, which
is one quadrature or one pair of cdf evaluations per component.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- y:

  A numeric vector of observations inside the truncation interval.
  Outside it the parent's value is returned, the method testing nothing.

- theta:

  A named list of parameters of the **parent**, aligned by the generic.
  Truncation adds no parameter, the endpoints being constants.

- scale:

  Handled by the generic after dispatch; this method always returns the
  parameter scale, and passes `"parameter"` down so the chain rule is
  applied once rather than twice.

- ...:

  Passed to the parent's method.

## Value

A named list with one numeric vector per parameter, keyed by the
parent's `params`, each of length `length(y)`.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrapper;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic;
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
which does pay for the normalizing constant.

## Examples

``` r
d <- gaussian1_distrib()
tr <- truncated(d, lower = -2, upper = 3)
y <- c(-1, 0, 2)
th <- list(mu = 0.3, sigma = 1.4)

# Identical to the parent's, component for component.
all.equal(distrib_cross_y(tr, y, th), distrib_cross_y(d, y, th))
#> [1] TRUE

# The gradient is not: there log Z does depend on theta.
all.equal(distrib_gradient(tr, y, th), distrib_gradient(d, y, th))
#> [1] "Component “mu”: Mean relative difference: 0.05596171"  
#> [2] "Component “sigma”: Mean relative difference: 0.5870033"
```
