# Variance of the Generalized Gamma Distribution

Closed form from the first two raw moments, \\\operatorname{Var}(Y) =
m_2 - m_1^2\\ with \\m_k = a^k\\\Gamma\\(d+k)/p\\/\Gamma(d/p)\\. The
scale enters as a square and the two shapes through the gamma ratios.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive, each a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length equal to the longest of the
three components.

## Details

The difference is not formed. With \\k = d/p\\ and \\h = 1/p\\, the
cumulant generating function of \\\log(Y/m_1)\\ is \\K(s) = \sum\_{n \ge
2} \psi^{(n-1)}(k)\\h^n (s^n - s)/n!\\, and with \\e^{K(s)} = \sum_n a_n
s^n\\ the central moments of \\Y/m_1\\ are \$\$\mu_r = r!\sum\_{n \ge r}
S(n, r)\\a_n,\$\$ \\S\\ the Stirling numbers of the second kind, a sum
in which the leading orders that cancel in \\m_2 - m_1^2\\ are absent.
It is evaluated in double-double arithmetic where \\8h \le k\\, and as
\\\exp(\log m_2 - 2\log m_1) - 1\\ from Stirling's form of the log-gamma
differences elsewhere. Then \\\operatorname{Var}(Y) = m_1^2\mu_2\\.

## Notation

\\a \> 0\\ is the scale, \\d \> 0\\ and \\p \> 0\\ the two shapes, and
\\m_k\\ the \\k\\-th raw moment.

## See also

[`mean.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.GenGamma1Distrib.md),
[`skewness.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GenGamma1Distrib.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## Examples

``` r
d <- gengamma1_distrib()

# At p = 1 the family is a gamma of shape d, whose variance is a^2 d.
all.equal(variance(d, list(a = 2, d = 3, p = 1)), 12)
#> [1] TRUE

# At d = p it agrees with the Weibull.
all.equal(variance(d, list(a = 2, d = 3, p = 3)),
          variance(weibull1_distrib(), list(mu = 2, sigma = 3)))
#> [1] TRUE

# Towards the lognormal the coefficient of variation is about sqrt(1/d).
variance(d, list(a = 1, d = 1e6, p = 1)) / mean(d, list(a = 1, d = 1e6, p = 1))^2
#> [1] 1e-06
```
