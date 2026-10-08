# Weibull Distribution-Function Derivatives, Mean and Shape

The derivatives of \\F(q) = 1 - \exp\\-(q/b)^\sigma\\\\ in \\(m,
\sigma)\\ to order four, each a closed form from its own compiled
kernel, carried to the upper tail and the log scale by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).

## Arguments

- distrib:

  A `Weibull3Distrib` object.

- q:

  A numeric vector of strictly positive quantiles.

- theta:

  A named list with components `mean` and `sigma`.

- lower.tail, log:

  As in
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md).

- ...:

  Unused.

## Value

A named list, one component per parameter or per multi-index.

## See also

[`distrib_cdf.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull3Distrib.md).

## Examples

``` r
d <- weibull3_distrib()
distrib_grad_cdf(d, c(1, 3), list(mean = 4, sigma = 1.7), log = FALSE)
#> $mean
#> [1] -0.03067729 -0.12954737
#> 
#> $sigma
#> [1] -0.1132162 -0.1432255
#> 
```
