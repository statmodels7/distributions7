# Folded Response Gradient

Computes the derivative of the folded log-density in the response,
\$\$\frac{\partial \ell}{\partial x} = w\\ g(x) - (1-w)\\ g(-x),\$\$
with \\g = \partial\log f/\partial y\\ the parent's own response
gradient and \\w\\ the weight of the positive preimage. The MINUS sign
is the chain rule through the reflected preimage: moving \\x\\ up moves
\\-x\\ down.

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

\\f\\ is the parent's density, \\g\\ its response gradient, \\w\\ the
weight of the positive preimage and \\\ell\\ the folded log-density.

## See also

[`distrib_hess_y.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.FoldedDistrib.md)
for the second order,
[`distrib_gradient.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.FoldedDistrib.md)
for the derivatives in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
y <- c(0.2, 1, 3)

distrib_grad_y(d, y, theta)
#> [1] -0.1148149 -0.5785036 -1.8130051

# Against a numerical derivative of the folded log-density.
vapply(y, function(v)
  numDeriv::grad(function(z) distrib_pdf(d, z, theta, log = TRUE), v),
  numeric(1))
#> [1] -0.1148149 -0.5785036 -1.8130051

# The minus sign written out: the reflected preimage pushes the other way.
g0 <- gaussian1_distrib()
w <- distributions7:::fold_parts(g0, y, theta)$w
all.equal(distrib_grad_y(d, y, theta),
          w * distrib_grad_y(g0, y, theta) -
            (1 - w) * distrib_grad_y(g0, -y, theta))
#> [1] TRUE
```
