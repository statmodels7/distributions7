# Folded Response Hessian

Computes the second derivative of the folded log-density in the
response, \$\$\frac{\partial^2 \ell}{\partial x^2} = w\left(h(x) +
g(x)^2\right) + (1-w)\left(h(-x) + g(-x)^2\right) - \left(w g(x) - (1-w)
g(-x)\right)^2,\$\$ with \\g\\ and \\h\\ the parent's first and second
response derivatives. The last term is the square of
[`distrib_grad_y.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.FoldedDistrib.md)'s
value, so the expression is the mixture's second-order ratio less the
square of its first: the same moment-to-cumulant step the parameter
derivatives take, applied in the response.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  A numeric vector of observations, non-negative.

- theta:

  A named list of the parent's parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Notation

\\f\\ is the parent's density, \\g\\ and \\h\\ its first and second
response derivatives, \\w\\ the weight of the positive preimage and
\\\ell\\ the folded log-density.

## See also

[`distrib_grad_y.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.FoldedDistrib.md)
for the first order,
[`distrib_hessian.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.FoldedDistrib.md)
for the parameter derivatives, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
y <- c(0.2, 1, 3)

distrib_hess_y(d, y, theta)
#> [1] -0.5744607 -0.5873235 -0.6469585

# Against a numerical second derivative.
vapply(y, function(v)
  numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
                    v)[1, 1],
  numeric(1))
#> [1] -0.5744607 -0.5873235 -0.6469585

# Far from zero the fold is invisible and the curvature is the parent's,
# which for a gaussian is -1 / sigma^2.
c(folded = distrib_hess_y(d, 8, theta), parent = -1 / 1.2^2)
#>     folded     parent 
#> -0.6925944 -0.6944444 
```
