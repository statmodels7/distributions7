# Binomial Distribution Object

Creates a distribution object for the Binomial distribution
parameterized by the probability of success \\\mu\\ and a number of
trials \\n\\ (size).

## Usage

``` r
binomial_distrib(link_mu = logit_link(), size = 1)
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\ (probability).
  Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html)
  to ensure the parameter stays within (0, 1).

- size:

  Integer or Numeric Vector. The number of trials \\n\\. Can be a single
  scalar (default 1) or a vector of the same length as the observations
  \\y\\.

## Value

An S7 object of class `BinomialDistrib` (inheriting from
`discrete_distrib`) representing the Binomial distribution.

## Details

The Binomial distribution models the number of successes in \\n\\
independent trials, each with success probability \\\mu\\. The number of
trials \\n\\ is fixed in the constructor (`size`) and treated as known.

**Probability mass function:** \$\$P(Y=y; \mu, n) = \dbinom{n}{y} \mu^y
(1-\mu)^{n-y}, \quad y \in \\0, 1, \dots, n\\\$\$

**Cumulative distribution function:** \$\$F(q; \mu, n) =
\sum\_{k=0}^{\lfloor q \rfloor} \dbinom{n}{k} \mu^k (1-\mu)^{n-k}\$\$

**Quantile function:** the generalized inverse \\Q(p) = \min\\y : F(y)
\ge p\\\\.

**Score, observed and expected Hessian** (with respect to \\\mu\\):
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - n\mu}{\mu(1-\mu)},
\qquad \dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} -
\dfrac{n-y}{(1-\mu)^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{n}{\mu(1-\mu)}\$\$

**Moments:** mean \\n\mu\\, variance \\n\mu(1-\mu)\\, skewness
\\(1-2\mu)/\sqrt{n\mu(1-\mu)}\\, excess kurtosis \\(1 -
6\mu(1-\mu))/(n\mu(1-\mu))\\.

**Parameter domains:**

- \\\mu \in (0, 1)\\

- \\n \in \mathbb{Z}^+\\ (fixed in the constructor, may vary per
  observation)

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
are also available.

## See also

- [`distrib_pdf.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BinomialDistrib.md)
  for the probability mass function.

- [`distrib_cdf.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BinomialDistrib.md)
  for the quantile function.

- [`distrib_rng.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BinomialDistrib.md)
  for random number generation.

- [`distrib_gradient.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BinomialDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BinomialDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md)
  for the analytical expected Hessian.
