# Mean of the Elastic-Net Distribution

Returns \\\mu\\, the location. The density is symmetric about it because
both factors of the product are symmetric about the same point. The
location is therefore the mean, the median and the mode at once.

The value is recycled to the length that the three parameters imply, so
the result has one mean per observation when the components of `theta`
vary by observation.

## Arguments

- x:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- theta:

  A named list with components `mu`, `lambda` and `alpha`. It is aligned
  and validated by name, so a missing or out-of-bounds component signals
  an error.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of the length that the recycled parameters
imply. The value is the location; the two rates do not enter it, so a
setting that varies one of them repeats one number.

## See also

[`variance.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.EnetDistrib.md)
and
[`skewness.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.EnetDistrib.md)
for the other two closed-form moments,
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the fourth, which is not closed form here, and
[`distrib_quantile.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.EnetDistrib.md),
which confirms the median.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 3, lambda = 2, alpha = 0.5)

# Mean, median and mode coincide.
c(mean = mean(d, th), median = distrib_quantile(d, 0.5, th))
#>   mean median 
#>      3      3 

# The result has one value per observation when a parameter varies.
mean(d, list(mu = c(0, 3, 7), lambda = 2, alpha = 0.5))
#> [1] 0 3 7
```
