# Mean of the Generalized Gamma Distribution

Closed form: \\E\[Y\] = a\\\Gamma\\(d+1)/p\\/\Gamma(d/p)\\, the first
raw moment. None of the three parameters is the mean: \\a\\ is a scale
and \\d\\ and \\p\\ are shapes, and the gamma ratio is what carries the
shapes into the mean.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive, each a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length equal to the longest of the three
components.

## Details

The family nests four the toolkit ships separately, and each is a check
on this formula: the gamma at \\p = 1\\, the Weibull at \\d = p\\, the
exponential at \\d = p = 1\\ and the half-normal at \\a = \sqrt2, d = 1,
p = 2\\.

## Notation

\\a \> 0\\ is the scale, \\d \> 0\\ and \\p \> 0\\ the two shapes.

## See also

[`variance.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.GenGamma1Distrib.md),
[`gengamma_raw_moments()`](https://statmodels7.github.io/distributions7/reference/gengamma_raw_moments.md)
for the shared quantities,
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

## Examples

``` r
d <- gengamma1_distrib()

# At p = 1 the family is a gamma of shape d, whose mean is a d.
all.equal(mean(d, list(a = 2, d = 3, p = 1)), 6)
#> [1] TRUE

# At d = p it is a Weibull of scale a and shape p.
all.equal(mean(d, list(a = 2, d = 3, p = 3)),
          mean(weibull1_distrib(), list(mu = 2, sigma = 3)))
#> [1] TRUE

# At a = sqrt(2), d = 1, p = 2 it is the half-normal, of mean sqrt(2/pi).
all.equal(mean(d, list(a = sqrt(2), d = 1, p = 2)), sqrt(2 / pi))
#> [1] TRUE
```
