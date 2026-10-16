# The Distribution Function of a von Mises by Its Bessel Series

Returns \\F(x)\\ on \\\[-\pi, \pi)\\ from the Fourier expansion of the
density, integrated term by term, in place of a quadrature.

## Usage

``` r
vm_cdf(y, mu, kappa)
```

## Arguments

- y:

  A numeric vector of quantiles. Below \\-\pi\\ the value is 0 and at or
  above \\\pi\\ it is 1, so that `-Inf` gives 0 and `Inf` gives 1. `NA`
  and `NaN` give `NA`.

- mu:

  A numeric vector of mean directions, recycled to the length of `y`.

- kappa:

  A numeric vector of concentrations, recycled to the length of `y`. A
  value at or below zero gives `NA` for an angle inside the support.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of the recycled
length of the inputs.

## Why a series

The density has no elementary antiderivative, and the base class
integrates it numerically, with one quadrature per observation. The
series converges rapidly. From \\e^{\kappa\cos\theta} = I_0(\kappa) +
2\sum\_{j\ge1} I_j(\kappa) \cos(j\theta)\\, integrating from \\-\pi\\ to
\\x\\, \$\$F(x) = \frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)}
\sum\_{j\ge1} \frac{I_j(\kappa)}{j} \big\[\sin(j(x - \mu)) +
\sin(j(\pi + \mu))\big\].\$\$ The second sine is the lower limit, and it
makes this the distribution function of the family as written: the
support is \\\[-\pi, \pi)\\ with the direction inside it, not a variable
wrapped around the circle.

Only the ratios \\I_j/I_0\\ are needed, and
[`numericals7::bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
gives them by a backward recurrence whose loop runs over the series
index and not over the data. The series is cheap for that reason: the
\\n\\ quadratures are replaced by a few dozen vectorized steps.

## Number of terms

The series is truncated at \\\max\\20, \lceil 8.5\sqrt{\kappa}\rceil +
10\\\\ terms, with \\\kappa\\ the largest concentration among the
observations evaluated, which exceeds the number of terms at which the
partial sums agree to double precision.

The sum is accumulated over blocks of observations, because the natural
expression forms an \\n \times m\\ matrix, which at a hundred thousand
points and a concentration of a hundred is already hundreds of
megabytes. The result is clamped to \\\[0, 1\]\\ to remove the rounding
error of the partial sums, whose absolute error is about \\10^{-16}\\;
[`vm_cdf_tail()`](https://statmodels7.github.io/distributions7/reference/vm_cdf_tail.md)
handles the tails below that resolution.

## See also

[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md),
which calls this;
[`numericals7::bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
for the recurrence; and
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
for the family.
