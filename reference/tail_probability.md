# The Probability of the Requested Tail

Internal helpers of
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).
`tail_cdf()` returns the probability of the requested tail at `q`, from
the family's own distribution function for that tail (so that the upper
tail is not formed as \\1 - F\\, which loses its digits once \\F\\ is
close to one), with the attribute `"tail"` set to `"lower"` or
`"upper"`. Called with the arguments that
[`log_tail_route()`](https://statmodels7.github.io/distributions7/reference/log_tail_route.md)
has stored in `log_tail_stash`, it returns the stored value.
`tail_probability()` returns the requested tail: the value itself when
it carries the attribute `"tail"` for that tail, and otherwise, the
value being read as the distribution function, \\F\\ or \\1 - F\\.

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

A numeric vector of probabilities; for `tail_cdf()`, with the attribute
`"tail"`.

## Examples

``` r
d <- logistic_distrib()
Fq <- distributions7:::tail_cdf(d, 40, list(mu = 0, sigma = 1), FALSE)
# The upper tail survives where 1 - F is zero.
c(from_family = distributions7:::tail_probability(Fq, FALSE),
  one_minus_F = 1 - plogis(40))
#>  from_family  one_minus_F 
#> 4.248354e-18 0.000000e+00 
```
