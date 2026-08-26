# Poisson Distribution

Builds the distribution object for the Poisson family parametrized by
its mean \\\mu \> 0\\, with mass \\e^{-\mu}\mu^y/y!\\ on the
non-negative integers. The returned object carries closed-form
derivatives of the log-mass to fourth order, observed and expected, and
closed-form moments.

The family is **equidispersed**: its variance equals its mean at every
parameter value. Real counts are usually more variable than that, so a
Poisson fit is often the null a count analysis argues against;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
add a dispersion parameter and both contain this family in a limit.

## Usage

``` r
poisson_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and is the **canonical** link
  of this family.

## Value

An S7 object of class `PoissonDistrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"poisson"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`,
`params_interpretation` `c(mu = "mean")`, `n_params` `1`,
`params_bounds` the list of \\(0, \infty)\\, and `link_params` the one
link given here.

## The parametrization

The mass on \\y \in \\0, 1, 2, \dots\\\\ is \$\$P(Y = y; \mu) =
\dfrac{e^{-\mu}\mu^{y}}{y!},\$\$ with \\\mu \in (0, \infty)\\. The mean
and the variance are both \\\mu\\, the skewness is \\\mu^{-1/2}\\ and
the excess kurtosis \\\mu^{-1}\\, so the distribution is right skewed at
a small mean and close to Gaussian at a large one.

## The canonical link, and what it buys

The log-mass is \\y\log\mu - \mu - \log y!\\, which is linear in \\y\\
once \\\eta = \log\mu\\ is the parameter. The logarithm is therefore the
canonical link, and on its scale \$\$\dfrac{\partial \ell}{\partial
\eta} = y - \mu, \qquad \dfrac{\partial^2 \ell}{\partial \eta^2} =
-\mu.\$\$ The second of these carries no data, so the observed and the
expected information coincide: Fisher scoring and Newton's method take
the same step, and the iteratively reweighted least squares of a Poisson
log-linear model is both at once. Under any other link the two differ
and the choice matters.

## Derivatives

On the parameter scale, with \\\ell = y\log\mu - \mu - \log y!\\,
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\mu}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}, \qquad
\dfrac{\partial^3 \ell}{\partial \mu^3} = \dfrac{2y}{\mu^3}, \qquad
\dfrac{\partial^4 \ell}{\partial \mu^4} = -\dfrac{6y}{\mu^4},\$\$ and
their expectations follow by substituting \\\mathbb{E}\[Y\] = \mu\\:
\\0\\, \\-1/\mu\\, \\2/\mu^2\\ and \\-6/\mu^3\\.

## Estimation

The maximum likelihood estimate is the sample mean, \\\hat\mu = \bar
y\\, in closed form. A quick check of the equidispersion assumption is
to compare the sample variance with it; the example below does that on
data drawn from the family and on data that is not.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \> 0\\ the mean,
which is also the variance. \\\eta = \log\mu\\ is the parameter on the
link scale. The **canonical** link of an exponential family is the one
that makes the log-mass linear in the sufficient statistic.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Chapter 4. Wiley, Hoboken.

## See also

[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for overdispersed counts, both of which contain this family in a limit;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
of which it is the many-trials limit;
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for excess zeros;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
for a memoryless count law;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the mean;
[PoissonDistrib](https://statmodels7.github.io/distributions7/reference/PoissonDistrib.md)
for the class.

## Examples

``` r
d <- poisson_distrib()
d
#> Distribution: Poisson
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: log        | Domain: (0, Inf)

# The mass is R's own and sums to one over the support.
y <- c(0, 2, 7)
all.equal(distrib_pdf(d, y, list(mu = 3)), dpois(y, 3))
#> [1] TRUE
sum(distrib_pdf(d, 0:200, list(mu = 3)))
#> [1] 1

# Equidispersed, and closer to symmetric as the mean grows.
vapply(c(0.5, 3, 20), function(m) {
  th <- list(mu = m)
  c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
}, numeric(3))
#>          [,1]      [,2]       [,3]
#> mean 0.500000 3.0000000 20.0000000
#> var  0.500000 3.0000000 20.0000000
#> skew 1.414214 0.5773503  0.2236068

# The estimate is the sample mean, in closed form.
set.seed(4)
z <- distrib_rng(d, 1000, list(mu = 4.2))
c(fitted = unname(coef(fit_distrib(d, z))), sample_mean = mean(z))
#>      fitted sample_mean 
#>       4.087       4.087 

# Equidispersion is a testable assumption. Drawn from the family the
# sample variance matches the mean; drawn with extra variability it does
# not, and a negative binomial is the family to reach for.
set.seed(6)
over <- rpois(1000, lambda = rgamma(1000, shape = 2, rate = 2 / 4.2))
rbind(poisson = c(mean = mean(z), var = var(z)),
      overdispersed = c(mean = mean(over), var = var(over)))
#>                mean       var
#> poisson       4.087  4.149581
#> overdispersed 4.060 11.866266
```
