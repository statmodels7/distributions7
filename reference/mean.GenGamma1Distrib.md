# Mean of the Generalized Gamma Distribution

The mean has the closed form \\E\[Y\] =
a\\\Gamma\\(d+1)/p\\/\Gamma(d/p)\\, the first raw moment. The mean is a
function of all three parameters: \\a\\ is a scale, \\d\\ and \\p\\ are
shapes, and the gamma ratio carries the shapes into the mean.

## Arguments

- x:

  A `GenGamma1Distrib`, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- theta:

  A named list with components `a`, `d` and `p`, all positive; each is a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, whose length is the one that the recycled
parameters imply.

## Details

The logarithm of the gamma ratio is formed in Stirling's form, \\(k -
1/2)\log(1 + h/k) + h\log(k + h) - h\\ plus the difference of the two
remainders, for \\k = d/p \ge 10\\ and \\h = 1/p\\, where the difference
of the two log-gamma values would lose digits.

The family nests the gamma at \\p = 1\\, the Weibull at \\d = p\\ and
the exponential at \\d = p = 1\\, which the toolkit ships as families,
and the half-normal at \\a = \sqrt2, d = 1, p = 2\\, which the toolkit
builds as `fixed(folded(gaussian1_distrib()), mu = 0)`. Each is a check
on this formula.

## Notation

\\a \> 0\\ is the scale, \\d \> 0\\ and \\p \> 0\\ the two shapes.

## See also

[`variance.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.GenGamma1Distrib.md),
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
