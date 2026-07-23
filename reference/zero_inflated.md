# Zero-Inflated Distribution Object (Discrete)

Creates a zero-inflated version of an existing **discrete** distribution
object: a mixture of a point mass at zero (with probability \\\zeta\\,
parameter `zi`) and the original count distribution.

## Usage

``` r
zero_inflated(distrib, link_zi = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` whose support includes 0.

- link_zi:

  A link function object for the zero-inflation probability \\\zeta\\.
  Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html).

## Value

An S7 object of class `ZeroInflatedDistrib` (inheriting from
`discrete_distrib`).

## Details

\$\$ P(Y=y; \theta, \zeta) = \begin{cases} \zeta + (1 - \zeta)f(0;
\theta) & y = 0 \\ (1 - \zeta)f(y; \theta) & y \> 0 \end{cases} \$\$

The resulting object supports the full `distrib` API: pdf, cdf,
quantile, rng, analytical gradient, observed and expected Hessian (all
derived from the parent's), plus numerical moments via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).

## Examples

``` r
if (FALSE) { # \dontrun{
zip <- zero_inflated(poisson_distrib())
distrib_pdf(zip, 0:5, list(mu = 3, zi = 0.2))
} # }
```
