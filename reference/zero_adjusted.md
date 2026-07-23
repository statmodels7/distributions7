# Zero-Adjusted Distribution Object

Creates a zero-adjusted version of an existing distribution object,
automatically dispatching on its type:

- **Discrete** (support including 0): a hurdle model, where the
  probability at zero is replaced by \\\pi\\ and positive values follow
  the zero-truncated parent.

- **Continuous**: a mixed distribution with a point mass \\\pi\\ at zero
  and the original continuous density scaled by \\1-\pi\\.

The mixture probability is exposed as parameter `za`.

## Usage

``` r
zero_adjusted(distrib, link_za = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` or `continuous_distrib`.

- link_za:

  A link function object for the zero probability \\\pi\\. Defaults to
  [`logit_link`](https://rdrr.io/pkg/linkfunctions7/man/logit_link.html).

## Value

An S7 object of class `ZeroAdjustedDiscreteDistrib` or
`ZeroAdjustedContinuousDistrib`.

## Examples

``` r
if (FALSE) { # \dontrun{
zap <- zero_adjusted(poisson_distrib())
distrib_pdf(zap, 0:5, list(mu = 3, za = 0.3))

zagamma <- zero_adjusted(gamma_distrib())
distrib_rng(zagamma, 10, list(mu = 2, sigma2 = 1, za = 0.3))
} # }
```
