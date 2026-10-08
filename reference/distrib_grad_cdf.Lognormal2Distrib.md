# Lognormal Distribution-Function Derivatives, Mean and Variance of Y

The derivatives of \\F(q) = \Phi\\(\log q - \mu_l)/\sqrt{S}\\\\ in \\(m,
v)\\ to order four, each a closed form from its own compiled kernel,
carried to the upper tail and the log scale by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).

## Arguments

- distrib:

  A `Lognormal2Distrib` object.

- q:

  A numeric vector of strictly positive quantiles.

- theta:

  A named list with components `mean` and `var`.

- lower.tail, log:

  As in
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md).

- ...:

  Unused.

## Value

A named list, one component per parameter or per multi-index.

## See also

[`distrib_cdf.Lognormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md).

## Examples

``` r
d <- lognormal2_distrib()
distrib_grad_cdf(d, c(1, 3), list(mean = 3, var = 2), log = FALSE)
#> $mean
#> [1] -0.05170637 -0.31582145
#> 
#> $var
#> [1] 0.02019204 0.01973884
#> 
```
