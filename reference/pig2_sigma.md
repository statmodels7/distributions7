# The Dispersion a Poisson-Inverse Gaussian Alpha Implies

Converts the orthogonal parametrization's \\\alpha\\ into the dispersion
\\\sigma\\ of
[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md):
\\\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2\\, the positive root
of \\\alpha^2\sigma^2 - 2\mu\sigma - 1 = 0\\, with no cancellation
anywhere in the domain.

## Usage

``` r
pig2_sigma(mu, alpha)
```

## Arguments

- mu:

  The mean, a positive numeric vector.

- alpha:

  The orthogonal parameter, a positive numeric vector.

## Value

A numeric vector of dispersions.

## See also

[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
