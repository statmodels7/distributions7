# Generalized Gamma Distribution Object

Builds a generalized gamma distribution object in Stacy's form, with a
scale \\a \> 0\\ and two shapes \\d \> 0\\ and \\p \> 0\\. It is the
flexible family for a positive response, and it makes the choice between
the gamma, the Weibull and the exponential something to estimate.

## Usage

``` r
gengamma1_distrib(
  link_a = log_link(),
  link_d = log_link(),
  link_p = log_link()
)
```

## Arguments

- link_a:

  A `linkfunctions7` link object for the scale \\a\\, which must be
  strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_d:

  A link object for the first shape \\d\\, also strictly positive.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_p:

  A link object for the second shape \\p\\, also strictly positive.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[GenGamma1Distrib](https://statmodels7.github.io/distributions7/reference/GenGamma1Distrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("a", "d", "p")`, its `bounds` `c(0, Inf)`, and its `link_params` the
three links given here.

## Density

\$\$f(y) = \dfrac{p}{a^{d}\\\Gamma(d/p)}\\y^{d-1}e^{-(y/a)^{p}}, \qquad
y \> 0.\$\$ The three parameters do three separate things: \\a\\ sets
the scale, \\d\\ the behavior at the origin, and \\p\\ the weight of the
upper tail.

## What it nests

Stacy's parametrization is chosen to make these visible:

- \\p = 1\\ is the
  [gamma](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
  with shape \\d\\ and scale \\a\\;

- \\d = p\\ is the
  [Weibull](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
  with shape \\p\\ and scale \\a\\;

- \\d = p = 1\\ is the
  [exponential](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md);

- \\a = \sqrt2\\, \\d = 1\\, \\p = 2\\ is the half-normal;

- \\p \to 0\\ with \\d/p\\ held large approaches the
  [lognormal](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

The first four are exact and are checked in the examples. The fifth is a
limit, reached at no admissible value.

## Score and information

With \\w = (y/a)^{p}\\, \\L = \log(y/a)\\ and \\k = d/p\\,
\$\$\dfrac{\partial\ell}{\partial a} = \dfrac{pw-d}{a}, \qquad
\dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
\dfrac{\partial\ell}{\partial p} = \dfrac{1}{p} +
\dfrac{d\psi(k)}{p^{2}} - wL,\$\$ and the expected information is
**closed form**. The reason is the one representation the whole family
rests on: \\u = (Y/a)^p\\ is Gamma with shape \\k\\ and unit rate, so
every expectation the Hessian needs is one of \\E\[u\] = k\\, \\E\[u\log
u\] = k\psi(k+1)\\ and \\E\[u(\log u)^2\] = k\\\psi(k+1)^2 +
\psi'(k+1)\\\\. The same representation gives the distribution function,
the quantile function and the generator.

## Moments

\\E\[Y^{r}\] = a^{r}\Gamma\\(d+r)/p\\/\Gamma(d/p)\\, finite for every
\\r \> -d\\. Unlike the generalized Pareto's, they exist at every
parameter value.

## Identification

The three parameters are **weakly identified together**: \\d\\ and \\p\\
enter the density largely through their ratio, and the profile
likelihood in that direction is flat. Measured at \\a = 2, d = 3, p =
1.5\\, the information's eigenvalues are 4.72, 0.138 and 0.00321, a
condition number of 1470, with the flat direction \\(0.79, -0.52,
0.33)\\. A fit of all three wants several hundred observations: on 200
the standard errors are 0.88, 0.59 and 0.44 against estimates of 2.20,
2.89 and 1.68, and on 2000 they are 0.31, 0.21 and 0.14. Holding one
with
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
is often the better model.

## Parameter domains

- \\a \in (0, \infty)\\

- \\d \in (0, \infty)\\

- \\p \in (0, \infty)\\

## References

Stacy, E. W. (1962). A generalization of the gamma distribution. *Annals
of Mathematical Statistics* 33, 1187-1192.

## See also

[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
and
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for the families it nests or approaches,
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
for holding a shape, and
[GenGamma1Distrib](https://statmodels7.github.io/distributions7/reference/GenGamma1Distrib.md)
for the class and its method list.

## Examples

``` r
d <- gengamma1_distrib()
d@params
#> [1] "a" "d" "p"
th <- list(a = 2, d = 3, p = 1.5)

distrib_pdf(d, c(0.5, 2, 5), th)
#> [1] 0.04136704 0.27590958 0.08999981

# The four exact special cases, at three points each.
y <- c(0.5, 2, 5)
c(gamma = max(abs(distrib_pdf(d, y, list(a = 2, d = 3, p = 1)) -
                  dgamma(y, shape = 3, scale = 2))),
  weibull = max(abs(distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)) -
                    dweibull(y, shape = 1.5, scale = 2))),
  exponential = max(abs(distrib_pdf(d, y, list(a = 2, d = 1, p = 1)) -
                        dexp(y, rate = 1 / 2))),
  half_normal = max(abs(distrib_pdf(d, y, list(a = sqrt(2), d = 1, p = 2)) -
                        2 * dnorm(y))))
#>        gamma      weibull  exponential  half_normal 
#> 2.775558e-17 5.551115e-17 2.775558e-17 1.110223e-16 

# The moment formula.
c(closed = 2 * gamma((3 + 1) / 1.5) / gamma(3 / 1.5), ours = mean(d, th))
#>   closed     ours 
#> 3.009151 3.009151 

# Two fits at two sizes: the standard errors say how much data the three
# parameters want.
for (n in c(200, 2000)) {
  set.seed(52)
  f <- fit_distrib(d, distrib_rng(d, n, th))
  print(rbind(estimate = coef(f), se = sqrt(diag(vcov(f)))))
}
#>                  a         d         p
#> estimate 2.2017942 2.8938782 1.6759546
#> se       0.8790896 0.5919983 0.4428119
#>                  a         d         p
#> estimate 2.0409360 3.0481613 1.5475668
#> se       0.3061717 0.2117766 0.1350711
```
