# The Concentration Above Which the Beta-Binomial Uses the Direct Sums

Returns 0.1, the dispersion \\\sigma\\ of
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
below which the derivatives are computed from the direct sums of
[`betabinom_direct()`](https://statmodels7.github.io/distributions7/reference/betabinom_direct.md)
rather than through the shapes, and whose reciprocal 10 is the total
shape \\\alpha + \beta\\ of
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
above which the third and fourth derivatives take the same route. The
value equals `bb1_direct_cut()` in the compiled kernels, which apply it
to the first and second derivatives. Through the shapes the polygamma
differences pass through map derivatives of order \\\sigma^{-(K+1)}\\
(betabinom1) or cancel among themselves (betabinom2) as the
concentration grows; above 0.1 that route is within \\10^{-13}\\ at
every order.

## Usage

``` r
bb_direct_cut()
```

## Value

A single number.

## See also

[`betabinom_direct()`](https://statmodels7.github.io/distributions7/reference/betabinom_direct.md).
