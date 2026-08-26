# Exponential Distribution

Builds the distribution object for the exponential family parametrized
by its **mean** \\\mu \> 0\\, with density \\\mu^{-1}e^{-y/\mu}\\ on
\\\[0, \infty)\\. The returned object carries closed-form derivatives of
the log-density to fourth order and closed-form moments.

The family has one parameter and a fixed shape: the standard deviation
equals the mean, the skewness is 2 and the excess kurtosis is 6 whatever
\\\mu\\ is. It is the only continuous law with a constant hazard, which
is the memorylessness the details set out.

## Usage

``` r
exponential_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `ExponentialDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"exponential"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
`params_interpretation` `c(mu = "mean")`, `n_params` `1`,
`params_bounds` the list of \\(0, \infty)\\, and `link_params` the one
link given here.

## The parametrization

The density on \\y \in \[0, \infty)\\ is \$\$f(y; \mu) =
\dfrac{1}{\mu}\exp\left(-\dfrac{y}{\mu}\right),\$\$ with \\\mu \in (0,
\infty)\\. The distribution function is \\F(q) = 1 - e^{-q/\mu}\\ and
the quantile function \\Q(p) = -\mu\log(1-p)\\.

The mean and the standard deviation are both \\\mu\\, so the coefficient
of variation is 1; the skewness is 2 and the excess kurtosis 6, neither
depending on the parameter. The median is \\\mu\log 2\\, below the mean.

R parametrizes its own `dexp` by the **rate**, \\1/\mu\\, and the
methods convert.
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
is the corresponding choice for the Laplace, where the rate is what the
package carries.

## Memorylessness

The survival function is \\P(Y \> q) = e^{-q/\mu}\\, so \\P(Y \> s + t
\mid Y \> s) = P(Y \> t)\\ for every \\s, t \ge 0\\: a wait already
endured tells nothing about the wait remaining. Equivalently the hazard
\\f/(1-F)\\ is the constant \\1/\mu\\. On the log scale this is the
statement that both \\\log f\\ and \\\log(1-F)\\ are straight lines in
\\y\\, which is why
[`distrib_grad_y.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ExponentialDistrib.md)
is a constant and
[`distrib_hess_y.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ExponentialDistrib.md)
is zero.

The exponential is the only continuous law with this property, so a
fitted exponential is also the statement that the process has no memory.
When it does,
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
and
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
both contain this family at a unit shape and let the hazard rise or
fall.

## Derivatives

Writing \\\ell = -\log\mu - y/\mu\\, the four orders are
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\mu^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{\mu - 2y}{\mu^3},
\qquad \dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{6y -
2\mu}{\mu^4}, \qquad \dfrac{\partial^4 \ell}{\partial \mu^4} =
\dfrac{6\mu - 24y}{\mu^5},\$\$ and their expectations follow by
substituting \\\mathbb{E}\[Y\] = \mu\\: \\0\\, \\-1/\mu^2\\, \\4/\mu^3\\
and \\-18/\mu^4\\. The information is \\1/\mu^2\\, so the asymptotic
standard error of \\\hat\mu\\ is \\\mu/\sqrt n\\.

## Estimation

The maximum likelihood estimate is the sample mean, \\\hat\mu = \bar
y\\, in closed form.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches it on the link scale; the example below checks it.

## Notation

\\\ell\\ is the log-density of one observation and \\\mu \> 0\\ the
mean, which is also the standard deviation and the reciprocal of the
rate. The **hazard** is \\f(y)/(1 - F(y))\\, the instantaneous failure
rate given survival to \\y\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 19. Wiley, New
York.

## See also

[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
and
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
which contain this family at a unit shape and let the hazard vary;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md),
its memoryless discrete counterpart;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md),
which contains it at a zero shape;
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which produces it from a larger family by holding a shape;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the mean;
[ExponentialDistrib](https://statmodels7.github.io/distributions7/reference/ExponentialDistrib.md)
for the class.

## Examples

``` r
d <- exponential_distrib()
d
#> Distribution: Exponential
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: log        | Domain: (0, Inf)

# The density is R's own at rate = 1/mu.
y <- c(0.3, 1.1, 4.0)
all.equal(distrib_pdf(d, y, list(mu = 2)), dexp(y, rate = 1 / 2))
#> [1] TRUE

# A fixed shape: the standard deviation is the mean, and the skewness and
# excess kurtosis do not move with mu.
vapply(c(0.5, 2, 10), function(m) {
  th <- list(mu = m)
  c(sd = std_dev(d, th), skew = skewness(d, th), kurt = kurtosis(d, th))
}, numeric(3))
#>      [,1] [,2] [,3]
#> sd    0.5    2   10
#> skew  2.0    2    2
#> kurt  6.0    6    6

# Memoryless: the chance of surviving one more unit is the same at every
# elapsed time.
vapply(c(0, 1, 5, 20), function(s)
  distrib_cdf(d, s + 1, list(mu = 2), lower.tail = FALSE) /
    distrib_cdf(d, s, list(mu = 2), lower.tail = FALSE), numeric(1))
#> [1] 0.6065307 0.6065307 0.6065307 0.6065307

# The estimate is the sample mean, in closed form.
set.seed(21)
z <- distrib_rng(d, 1000, list(mu = 3))
c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#>      fitted sample_mean 
#>    2.821234    2.821234 

# It is a Weibull of unit shape, which fixed() produces from the larger
# family; the two densities agree exactly.
all.equal(distrib_pdf(d, y, list(mu = 2)),
          distrib_pdf(fixed(weibull1_distrib(), sigma = 1), y,
                      list(mu = 2)))
#> [1] TRUE
```
