# Skew Normal Random Generation

Draws from the skew normal exactly, using Azzalini's stochastic
representation: with \\U_0\\ and \\U_1\\ independent standard Gaussians
and \\\delta = \alpha/\sqrt{1+\alpha^2}\\, \$\$Z = \delta\\\|U_0\| +
\sqrt{1-\delta^2}\\U_1\$\$ is standard skew normal with shape
\\\alpha\\, and \\Y = \mu + \sigma Z\\. No inversion and no rejection is
involved. The cost is two `rnorm` calls whatever the shape is, and the
draws follow the density exactly.

The representation also reads off why the family cannot be very skewed:
the half-normal component enters with weight \\\delta\\, which saturates
at one as \\\alpha \to \infty\\, so the most skewed member is the
half-normal itself.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of length `n`; a component of length 1 is
  recycled, so a parameter varying by observation draws one value per
  observation from its own member of the family.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## Notation

\\\delta = \alpha/\sqrt{1+\alpha^2}\\ is the correlation between \\Z\\
and the half-normal component; \\\mu\\ is the location and \\\sigma\\
the scale.

## See also

[`distrib_pdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal1Distrib.md)
for the density the draws follow,
[`mean.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md)
and
[`variance.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md)
for the moments they should reproduce, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
th <- list(mu = 0, sigma = 1, alpha = 3)

set.seed(1)
x <- distrib_rng(d, 2e5, th)

# Three moments against their closed forms.
rbind(sample = c(mean(x), var(x), mean((x - mean(x))^3) / sd(x)^3),
      theory = c(mean(d, th), variance(d, th), skewness(d, th)))
#>             [,1]      [,2]      [,3]
#> sample 0.7583858 0.4279862 0.6619280
#> theory 0.7569398 0.4270422 0.6670236

# The representation, written out at the same seed.
delta <- 3 / sqrt(1 + 9)
set.seed(1)
z <- delta * abs(rnorm(2e5)) + sqrt(1 - delta^2) * rnorm(2e5)
all.equal(x, z)
#> [1] TRUE
```
