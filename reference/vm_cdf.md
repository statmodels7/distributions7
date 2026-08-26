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
  above \\\pi\\ it is 1. `NA` and non-finite values give `NA`.

- mu:

  A numeric vector of mean directions, recycled to the length of `y`.

- kappa:

  A numeric vector of concentrations, recycled to the length of `y`. A
  value at or below zero gives `NA`.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of the recycled
length of the inputs.

## Why a series

The density has no elementary antiderivative, and the base class
integrates it numerically: one quadrature per observation, which
measured 4.5 seconds at ten thousand points and made this family a
thousand times dearer than any other. The series is rapidly convergent.
From \\e^{\kappa\cos\theta} = I_0(\kappa) + 2\sum\_{j\ge1} I_j(\kappa)
\cos(j\theta)\\, integrating from \\-\pi\\ to \\x\\, \$\$F(x) =
\frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)} \sum\_{j\ge1}
\frac{I_j(\kappa)}{j} \big\[\sin(j(x - \mu)) + \sin(j(\pi +
\mu))\big\].\$\$ The second sine is the lower limit, and it makes this
the distribution function of the family **as written**: the support is
\\\[-\pi, \pi)\\ with the direction inside it, not a variable wrapped
around the circle.

Only the **ratios** \\I_j/I_0\\ are needed, and
[`numericals7::bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
gives them by a backward recurrence whose loop runs over the series
index and not over the data. The series is cheap for that reason: \\n\\
quadratures become a few dozen vectorized steps.

## How many terms, measured

Compared against the same series at four times the length, machine
precision is reached at 10 terms at \\\kappa = 0.5\\, 26 at 10, 90 at
100, 242 at 1000 and 404 at 3000, always under \\8.5\sqrt{\kappa} +
10\\, which is the rule used, with a floor of 20.

The sum is accumulated over blocks of observations, because the natural
expression forms an \\n \times m\\ matrix, which at a hundred thousand
points and a concentration of a hundred is already hundreds of
megabytes. The result is clamped to \\\[0, 1\]\\: the series is exact
and its rounding is not, and a distribution function cannot leave its
own interval.

## See also

[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md),
which calls this;
[`numericals7::bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
for the recurrence; and
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
for the family.
