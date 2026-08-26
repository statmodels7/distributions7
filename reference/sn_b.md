# The Constant Behind the Centered Parametrization

Returns \\b = \sqrt{2/\pi} \approx 0.7978846\\, which is \\E\[\|Z\|\]\\
for a standard Gaussian \\Z\\. It is the one constant the skew normal's
first moment introduces, and it appears in every quantity of the
centered parametrization: the mean is \\\xi + \omega b\delta\\, the
variance \\\omega^2(1-b^2\delta^2)\\, and the largest reachable skewness
is built from \\b/\sqrt{1-b^2}\\.

## Usage

``` r
sn_b()
```

## Value

A single number.

## See also

[`sn_max_skew()`](https://statmodels7.github.io/distributions7/reference/sn_max_skew.md),
which is written in terms of it, and
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the parametrization.

## Examples

``` r
distributions7:::sn_b()
#> [1] 0.7978846

# It is the mean of the absolute value of a standard Gaussian.
set.seed(1)
c(constant = distributions7:::sn_b(), sample = mean(abs(rnorm(1e6))))
#>  constant    sample 
#> 0.7978846 0.7981197 
```
