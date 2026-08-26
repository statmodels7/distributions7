# Geometric Distribution

Builds the distribution object for the geometric family parametrized by
its **mean** \\\mu \> 0\\: the number of failures before the first
success in independent trials of success probability \\p = 1/(1+\mu)\\.
The returned object carries closed-form derivatives of the log-mass to
fourth order, observed and expected, and closed-form moments.

The variance is \\\mu(1+\mu)\\, always above the mean, so a geometric
fit is an overdispersed alternative to a Poisson with no extra parameter
to estimate. It is the negative binomial at a dispersion of exactly one,
and it is the only discrete law that is memoryless.

## Usage

``` r
geometric_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `GeometricDistrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"geometric"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
`params_interpretation` `c(mu = "mean")`, `n_params` `1`,
`params_bounds` the list of \\(0, \infty)\\, and `link_params` the one
link given here.

## The parametrization

The mass on \\y \in \\0, 1, 2, \dots\\\\ is \$\$P(Y = y; \mu) =
\dfrac{1}{1+\mu}\left(\dfrac{\mu}{1+\mu}\right)^{y},\$\$ with \\\mu \in
(0, \infty)\\. The family is parametrized by the mean and not by the
success probability, so that a link can carry it to the unconstrained
scale;
[`geom_prob()`](https://statmodels7.github.io/distributions7/reference/geom_prob.md)
converts, giving \\p = 1/(1+\mu)\\ for the base R functions.

The mean is \\\mu\\, the variance \\\mu(1+\mu)\\, and the mode is always
0: consecutive masses differ by the constant factor \\\mu/(1+\mu)\\,
whatever the mean is. The distribution is right skewed at every
parameter value.

## Overdispersion and memorylessness

The variance exceeds the mean by \\\mu^2\\, so this family sits above a
Poisson of the same mean and is a one-parameter answer to
overdispersion. The price is precision: the information is
\\1/(\mu(1+\mu))\\ against a Poisson's \\1/\mu\\. Where the amount of
overdispersion should itself be estimated,
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
adds the parameter and contains this family at \\\theta = 1\\.

The survival function is \\P(Y \> q) = (\mu/(1+\mu))^{q+1}\\, so \\P(Y
\> s + t \mid Y \> s) = P(Y \> t)\\: the count already accumulated tells
nothing about the count remaining. The geometric is the only discrete
law with this property, as the exponential is the only continuous one,
and the two are counterparts.

## Estimation

The maximum likelihood estimate is the sample mean, \\\hat\mu = \bar
y\\, in closed form; the score is the residual divided by the variance
and its sum vanishes there.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \> 0\\ the mean,
with variance \\\mu(1+\mu)\\. The success probability is \\p =
1/(1+\mu)\\. Counting starts at zero: \\Y\\ is the number of failures
**before** the first success, not the number of trials.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Chapter 5. Wiley, Hoboken.

## See also

[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md),
of which this is the case \\\theta = 1\\;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the equidispersed alternative;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md),
the memoryless continuous counterpart;
[`geom_prob()`](https://statmodels7.github.io/distributions7/reference/geom_prob.md)
for the conversion to the success probability;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the mean;
[GeometricDistrib](https://statmodels7.github.io/distributions7/reference/GeometricDistrib.md)
for the class.

## Examples

``` r
d <- geometric_distrib()
d
#> Distribution: Geometric
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: log        | Domain: (0, Inf)

# The mass is R's own at prob = 1/(1+mu), and sums to one.
all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
          dgeom(c(0, 2, 7), prob = 1 / 4))
#> [1] TRUE
sum(distrib_pdf(d, 0:400, list(mu = 3)))
#> [1] 1

# Overdispersed by construction: the variance is mu(1+mu).
vapply(c(0.5, 3, 20), function(m) {
  th <- list(mu = m)
  c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
}, numeric(3))
#>          [,1]      [,2]       [,3]
#> mean 0.500000  3.000000  20.000000
#> var  0.750000 12.000000 420.000000
#> skew 2.309401  2.020726   2.000595

# Memoryless, as the exponential is on the continuous side.
vapply(c(0, 2, 10, 50), function(s)
  distrib_cdf(d, s + 1, list(mu = 3), lower.tail = FALSE) /
    distrib_cdf(d, s, list(mu = 3), lower.tail = FALSE), numeric(1))
#> [1] 0.75 0.75 0.75 0.75

# It is a negative binomial at a dispersion of one.
all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
          distrib_pdf(fixed(negbin2_distrib(), theta = 1), c(0, 2, 7),
                      list(mu = 3)))
#> [1] TRUE

# The estimate is the sample mean, in closed form.
set.seed(5)
z <- distrib_rng(d, 2000, list(mu = 2.5))
c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#>      fitted sample_mean 
#>      2.4455      2.4455 
```
