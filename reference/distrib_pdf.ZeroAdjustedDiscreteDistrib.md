# Zero-Adjusted Discrete Probability Mass Function

Computes the hurdle mass function \$\$P(Y = 0) = \pi, \qquad P(Y = y) =
(1-\pi)\frac{f(y;\theta)}{1 - f(0;\theta)} \quad (y \> 0).\$\$ The mass
at zero is the parameter itself, so it may sit ABOVE or BELOW the
parent's \\f(0)\\; the positive part is the parent renormalized away
from zero. Everything is computed on the log scale and exponentiated at
the end, with the normalizing constant taken through
[`base::log1p()`](https://rdrr.io/r/base/Log.html) so that a parent with
almost no mass at zero does not lose digits to `log(1 - f0)`.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `za`, each a
  numeric vector of length 1 or of the length of `y`. `za` must lie
  strictly inside \\(0, 1)\\.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned, which is
  the quantity actually computed. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Notation

\\f\\ is the parent's mass function, \\F\\ its distribution function,
\\\pi\\ the probability of a zero, \\s\\ the parent's score and \\\ell\\
the log-mass of one observation.

## See also

[`distrib_cdf.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedDiscreteDistrib.md)
for the distribution function,
[`distrib_pdf.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroInflatedDistrib.md)
for the wrapper that adds to the mass at zero where this one replaces
it, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)

distrib_pdf(d, 0:5, theta)
#> [1] 0.40000000 0.09431225 0.14146838 0.14146838 0.10610129 0.06366077

# The mass at zero is the parameter, and the rest is the parent
# renormalized away from it.
c(at_zero = distrib_pdf(d, 0, theta), parameter = 0.4)
#>   at_zero parameter 
#>       0.4       0.4 
all.equal(distrib_pdf(d, 1:5, theta),
          0.6 * dpois(1:5, 3) / (1 - dpois(0, 3)))
#> [1] TRUE

# It can produce FEWER zeros than the parent, which inflation cannot.
c(hurdle = distrib_pdf(d, 0, list(mu = 3, za = 0.01)),
  parent = dpois(0, 3))
#>     hurdle     parent 
#> 0.01000000 0.04978707 

# And it is a distribution: the mass sums to one.
sum(distrib_pdf(d, 0:200, theta))
#> [1] 1
```
