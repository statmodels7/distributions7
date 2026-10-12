# The Probability of the Requested Tail

Internal helpers of
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).
`tail_cdf()` returns the distribution function at `q` and, when the
upper tail is wanted, attaches the family's own upper-tail probability
as the attribute `"upper"`, so that the upper tail is not formed as
\\1 - F\\, which loses its digits once \\F\\ is close to one.
`tail_probability()` returns the lower tail, or the upper tail from that
attribute when it is present and from \\1 - F\\ otherwise.

## Usage

``` r
tail_cdf(distrib, q, theta, lower.tail)

tail_probability(Fq, lower.tail)
```

## Arguments

- distrib:

  A distribution object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameter values.

- lower.tail:

  Is the lower tail wanted? A single logical.

- Fq:

  The value of `tail_cdf()`, or a plain numeric vector of distribution
  function values.

## Value

A numeric vector of probabilities; for `tail_cdf()`, the distribution
function with the attribute `"upper"` when `lower.tail` is `FALSE`.

## Examples

``` r
d <- logistic_distrib()
Fq <- distributions7:::tail_cdf(d, 40, list(mu = 0, sigma = 1), FALSE)
# The upper tail survives where 1 - F is zero.
c(from_family = distributions7:::tail_probability(Fq, FALSE),
  one_minus_F = 1 - as.numeric(Fq))
#>  from_family  one_minus_F 
#> 4.248354e-18 0.000000e+00 
```
