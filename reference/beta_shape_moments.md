# Moments of the Beta in Its Mean and Precision

Returns the variance, the skewness or the excess kurtosis of the beta
law with mean \\\mu\\ and precision \\\phi = \alpha + \beta\\, in forms
in which the shapes are divided out of every quotient:
\$\$\mathrm{Var}(Y) = \frac{\mu(1-\mu)}{1+\phi},\qquad \gamma_1 =
\frac{2(1-2\mu)\sqrt{1+\phi}}{(2+\phi)\sqrt{\mu(1-\mu)}},\$\$
\$\$\gamma_2 = \frac{6\[(1-2\mu)^2(1+\phi)/(2+\phi) - \mu(1-\mu)\]}
{\mu(1-\mu)(3+\phi)}.\$\$ Written in the shapes, the same quotients
carry \\\alpha\beta\\, which underflows below \\\phi \method 10^{-154}\\
and overflows above \\10^{154}\\; here \\\phi\\ enters only through
\\1 + \phi\\, \\2 + \phi\\ and \\3 + \phi\\, and never in a product of
two of them.

## Usage

``` r
beta_shape_moments(mu, phi, what)
```

## Arguments

- mu:

  The mean, a numeric vector in \\(0, 1)\\.

- phi:

  The precision, a positive numeric vector.

- what:

  One of `"variance"`, `"skewness"` or `"kurtosis"`.

## Value

A numeric vector of the recycled length of `mu` and `phi`.

## See also

[`skewness.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta1Distrib.md),
[`kurtosis.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Beta1Distrib.md)
and the moments of
[Beta2Distrib](https://statmodels7.github.io/distributions7/reference/Beta2Distrib.md),
which call it.

## Examples

``` r
distributions7:::beta_shape_moments(0.3, 1e-200, "kurtosis")
#> [1] -1.238095
```
