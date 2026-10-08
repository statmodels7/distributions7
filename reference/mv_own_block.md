# The Quantities a Structured Parametrization Is About

Returns the block
[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md)
declares, with the structure's common variance `scale` replaced by its
square root where the matrix is a covariance or a scale matrix, or
`NULL` where the structure declares no block. A structured matrix such
as `compound_symmetry()` or `ar1()` is fixed by these few quantities, so
they are what a summary reports in place of the standard deviation of
every coordinate and the correlation of every pair, which repeat them.

## Usage

``` r
mv_own_block(distrib, theta, sd_label = "sd")
```

## Arguments

- distrib:

  A multivariate distribution object.

- theta:

  Its parameters, aligned.

- sd_label:

  The name given to the square root of `scale`.

## Value

The result of
[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md),
possibly with one row renamed and transformed, or `NULL`.

## Details

The standard deviation is \\\sqrt{s}\\ with \\s\\ the `scale`, so its
Jacobian row is the `scale` row times \\1/(2\sqrt{s})\\; it is
intervalled on the log scale, as `scale` was. Where the matrix is a
precision the `scale` is a precision's common diagonal and is reported
as it is.

## See also

[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md)
for the declared block,
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
for the consumer.

## Examples

``` r
d <- mvgaussian1_distrib(3, parameters7::compound_symmetry(3))
th <- distributions7:::align_theta(d, list(mu1 = 0, mu2 = 0, mu3 = 0,
  sigma_log_scale = log(4), sigma_logit_rho = 0))
distributions7:::mv_own_block(d, th)$value
#>   sd  rho 
#> 2.00 0.25 
```
