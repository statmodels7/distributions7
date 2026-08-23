# The Distribution Function of a von Mises by Its Bessel Series

\\F(x)\\ on \\\[-\pi, \pi)\\, from the Fourier expansion of the density
integrated term by term.

## Usage

``` r
vm_cdf(y, mu, kappa)
```

## Arguments

- y:

  The quantiles.

- mu:

  The mean directions.

- kappa:

  The concentrations.

## Value

The distribution function at `y`.

## Details

The density has no elementary antiderivative, and the base class
integrates it numerically: one quadrature per observation, which
measured 4.5 seconds at ten thousand points and made this family a
thousand times dearer than any other. It has a rapidly convergent series
instead. From \\e^{\kappa\cos\theta} = I_0(\kappa) + 2\sum\_{j\ge1}
I_j(\kappa) \cos(j\theta)\\, integrating from \\-\pi\\ to \\x\\,
\$\$F(x) = \frac{x + \pi}{2\pi} + \frac{1}{\pi I_0(\kappa)}
\sum\_{j\ge1} \frac{I_j(\kappa)}{j} \big\[\sin(j(x - \mu)) +
\sin(j(\pi + \mu))\big\].\$\$ The second sine is the lower limit and is
what makes this the distribution function of the family as written – the
support is \\\[-\pi, \pi)\\ with the location inside it, not a variable
wrapped around the circle.

Only the RATIOS \\I_j/I_0\\ are needed and
[`bessel_i_ratios`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
gives them by a backward recurrence whose loop runs over the series
index rather than over the data. That is what makes the series cheaper:
\\n\\ quadratures become a few dozen vectorized steps.

HOW MANY TERMS is measured rather than assumed. Comparing against the
same series at four times the length, machine precision is reached at 10
terms at \\\kappa = 0.5\\, 26 at 10, 90 at 100, 242 at 1000 and 404 at
3000 – always under \\8.5\sqrt{\kappa} + 10\\, which is the rule used.
The sum is accumulated over blocks of observations because the natural
expression forms an \\n \times m\\ matrix, which at a hundred thousand
points and a concentration of a hundred is already hundreds of
megabytes.

## See also

[`bessel_i_ratios`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.html)
