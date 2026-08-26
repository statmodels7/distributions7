# Folded Density

Computes the folded density \$\$L(x; \theta) = f(x; \theta) + f(-x;
\theta), \qquad x \ge 0,\$\$ the parent's density at the two preimages
of \\x\\ under the absolute value, added together. Below zero the folded
variable places no mass, so the density is exactly `0` and its logarithm
`-Inf`.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  A numeric vector of observations. Negative values return `0` without
  an error, the folded support being \\\[0, \infty)\\.

- theta:

  A named list of the PARENT's parameters, each a numeric vector of
  length 1 or of the length of `y`. Folding adds and removes nothing.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. The
  logarithm is taken of the sum, not inside the parent, so it underflows
  where the sum does. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Details

At \\x = 0\\ the two preimages coincide and the folded density is twice
the parent's. The addition is what separates a fold from a change of
variable: the map is two to one, so there is no Jacobian to divide by,
and
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
cannot express it.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one and \\\theta\\ the
parameters shared by both.

## See also

[`distrib_cdf.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.FoldedDistrib.md)
for the distribution function,
[`distrib_gradient.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.FoldedDistrib.md)
for the score,
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for the family, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)

distrib_pdf(d, c(0, 0.5, 2), theta)
#> [1] 0.6096206 0.5673785 0.1901609

# Which is the parent's density at the two preimages, added.
y <- c(0, 0.5, 2)
all.equal(distrib_pdf(d, y, theta),
          dnorm(y, 0.5, 1.2) + dnorm(-y, 0.5, 1.2))
#> [1] TRUE

# At zero the preimages coincide, so the density is twice the parent's.
c(folded = distrib_pdf(d, 0, theta), twice = 2 * dnorm(0, 0.5, 1.2))
#>    folded     twice 
#> 0.6096206 0.6096206 

# Below zero there is no mass.
distrib_pdf(d, c(-1, -0.1), theta)
#> [1] 0 0

# And it integrates to one over the half line.
integrate(function(z) distrib_pdf(d, z, theta), 0, Inf)$value
#> [1] 1
```
