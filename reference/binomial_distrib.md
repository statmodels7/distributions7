# Binomial Distribution

Builds the distribution object for the binomial family: the number of
successes in `size` independent trials, each succeeding with probability
\\\mu \in (0, 1)\\. The returned object carries closed-form derivatives
of the log-mass to fourth order, observed and expected, and closed-form
moments.

**`size` is fixed data, not a parameter.** It is carried on the object,
it is not estimated, it has no link and no bound, and it does not appear
in `params`. Giving it one value per observation is how grouped binary
data with unequal group sizes is described.

The default link is the logit, the canonical link here, so the observed
and the expected information coincide on its scale.

## Usage

``` r
binomial_distrib(link_mu = logit_link(), size = 1)
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the probability \\\mu\\.
  Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  the canonical link, which maps \\(0, 1)\\ onto the whole line.
  [`linkfunctions7::probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.html)
  and
  [`linkfunctions7::cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.html)
  are the usual alternatives; under either the observed and expected
  information differ.

- size:

  A numeric vector of trial counts, positive integers. A single value
  applies to every observation; a vector of the length of the response
  gives one count per observation. Defaults to `1`, which makes the
  object a Bernoulli; pass the group sizes for grouped data. The value
  is stored in the `size` property and sets `bounds` to
  `c(0, max(size))`.

## Value

An S7 object of class `BinomialDistrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"binomial"`, `dimension`
`"univariate"`, `bounds` `c(0, max(size))`, `size` as given, `params`
`"mu"`, `params_interpretation` `c(mu = "probability")`, `n_params` `1`,
`params_bounds` the list of \\(0, 1)\\, and `link_params` the one link.

## The parametrization

The mass on \\y \in \\0, 1, \dots, n\\\\ is \$\$P(Y = y; \mu) =
\binom{n}{y}\mu^{y}(1-\mu)^{n-y},\$\$ with \\n\\ the number of trials
and \\\mu \in (0, 1)\\. The mean is \\n\mu\\, the variance
\\n\mu(1-\mu)\\, the skewness \\(1-2\mu)/\sqrt{n\mu(1-\mu)}\\ and the
excess kurtosis \\(1-6\mu(1-\mu))/(n\mu(1-\mu))\\. Both shape measures
shrink like \\1/n\\, which is the central limit theorem visible in the
moments.

As \\n \to \infty\\ with \\n\mu\\ held, the family tends to a Poisson of
that mean;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
is the limit and is the family to use when the number of trials is large
and unrecorded.

## size is data

The number of trials sits on the object rather than in `theta`, so it
takes no link, is not counted in `n_params`, and cannot be estimated.
Two consequences worth knowing: the same object cannot be reused across
data sets with different group sizes, and `bounds` records
`c(0, max(size))`, so a vector `size` gives a bound that is correct for
the largest group and loose for the rest.

## The canonical link

The log-mass is linear in \\y\\ once \\\eta = \log(\mu/(1-\mu))\\ is the
parameter, so the logit is canonical and on its scale
\$\$\dfrac{\partial \ell}{\partial \eta} = y - n\mu, \qquad
\dfrac{\partial^2 \ell}{\partial \eta^2} = -n\mu(1-\mu),\$\$ the second
carrying no data. Under a probit or a complementary log-log link the
observed and expected information differ.

## Estimation

The maximum likelihood estimate is the pooled proportion, \\\hat\mu =
\sum y_i / \sum n_i\\, in closed form. If the observed counts are more
variable than \\n\mu(1-\mu)\\ allows, the binomial is the wrong family
and
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
adds the dispersion parameter it lacks.

## Notation

\\\ell\\ is the log-mass of one observation, \\n\\ the number of trials
and \\\mu \in (0,1)\\ the success probability. \\\eta =
\log(\mu/(1-\mu))\\ is the log odds. The **canonical** link of an
exponential family makes the log-mass linear in the sufficient
statistic.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Chapter 3. Wiley, Hoboken.

## See also

[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md),
the case `size = 1`;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for overdispersed grouped binary data;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md),
the many-trials limit;
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
for more than two categories;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probability;
[BinomialDistrib](https://statmodels7.github.io/distributions7/reference/BinomialDistrib.md)
for the class.

## Examples

``` r
d <- binomial_distrib(size = 10)
d
#> Distribution: Binomial
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (probability)        | Link: logit      | Domain: (0, 1)

# The mass is R's own and sums to one over 0:n.
all.equal(distrib_pdf(d, c(0, 4, 10), list(mu = 0.3)),
          dbinom(c(0, 4, 10), size = 10, prob = 0.3))
#> [1] TRUE
sum(distrib_pdf(d, 0:10, list(mu = 0.3)))
#> [1] 1

# size = 1 is the Bernoulli, exactly.
all.equal(distrib_pdf(binomial_distrib(size = 1), c(0, 1), list(mu = 0.3)),
          distrib_pdf(bernoulli_distrib(), c(0, 1), list(mu = 0.3)))
#> [1] TRUE

# Skewness shrinks like 1/sqrt(n): the central limit theorem in the moments.
vapply(c(1, 10, 100), function(n)
  skewness(binomial_distrib(size = n), list(mu = 0.3)), numeric(1))
#> [1] 0.87287156 0.27602622 0.08728716

# Many trials at a small probability is a Poisson of the same mean.
rbind(binomial = distrib_pdf(binomial_distrib(size = 1000),
                             0:5, list(mu = 3 / 1000)),
      poisson = distrib_pdf(poisson_distrib(), 0:5, list(mu = 3)))
#>                [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> binomial 0.04956308 0.1491367 0.2241537 0.2243786 0.1682839 0.1008691
#> poisson  0.04978707 0.1493612 0.2240418 0.2240418 0.1680314 0.1008188

# Grouped data with unequal group sizes: one size per observation.
set.seed(5)
g <- binomial_distrib(size = c(5, 10, 20, 50))
rbind(size = c(5, 10, 20, 50), draw = distrib_rng(g, 4, list(mu = 0.3)))
#>      [,1] [,2] [,3] [,4]
#> size    5   10   20   50
#> draw    1    4    9   13
```
