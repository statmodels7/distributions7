# Poisson-Inverse Gaussian Probability Mass Function

Computes the Poisson-inverse Gaussian mass. With \\c = 1 + 2\sigma\mu\\
and \\\alpha = \sqrt{c}/\sigma\\, \$\$P(Y = y) =
\sqrt{\dfrac{2\alpha}{\pi}}\\ \dfrac{\mu^y
e^{1/\sigma}}{(\alpha\sigma)^y\\ y!}\\ K\_{y-1/2}(\alpha).\$\$

The value is not computed that way. At half-integer order the Bessel
function is a finite sum, the prefactors cancel, and what the compiled
kernel evaluates is \\\ell(y) = y\log\mu - (y/2)\log c + (1-\sqrt
c)/\sigma + \log S_y(\alpha) - \log y!\\ with \\S_y\\ a sum of \\y\\
positive terms taken on the log scale. Nothing cancels and no Bessel
routine is called.

## Arguments

- distrib:

  A `Pig1Distrib` object, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- y:

  A numeric vector of counts. A negative, non-integer or non-finite
  value is off the support and gives a probability of 0, or `-Inf` with
  `log = TRUE`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-probability is returned.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A numeric vector of probabilities, of the length of the recycled inputs.

## Notation

\\\mu\\ is the mean, \\\sigma\\ the dispersion, \\K\_\nu\\ the modified
Bessel function of the second kind, and \\S_y\\ the finite sum defined
in `pig1_pdf_cpp`.

## See also

[`distrib_gradient.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md)
for the score, `pig1_pdf_cpp` for the kernel,
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the same law in orthogonal coordinates, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
y <- 0:6
th <- list(mu = 3, sigma = 0.8)

# Against the Bessel formula written out, which the kernel does not use.
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
all.equal(distrib_pdf(d, y, th),
          sqrt(2 * al / pi) * 3^y * exp(1 / 0.8) /
            ((al * 0.8)^y * factorial(y)) * besselK(al, y - 0.5))
#> [1] TRUE

# The mass sums to one.
sum(distrib_pdf(d, 0:300, th))
#> [1] 1

# The upper tail is heavier than a negative binomial's at the same
# variance, which is the reason to prefer this family.
nb <- negbin2_distrib()
rbind(pig = distrib_pdf(d, c(20, 40, 60), th),
      negbin = distrib_pdf(nb, c(20, 40, 60), list(mu = 3, theta = 1 / 0.8)))
#>                [,1]         [,2]         [,3]
#> pig    0.0006790258 5.753016e-06 7.233866e-08
#> negbin 0.0004803803 5.368144e-07 5.596723e-10

# Off the support.
distrib_pdf(d, c(-1, 1.5), th)
#> [1] 0 0
```
