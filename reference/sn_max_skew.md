# The Largest Skewness a Skew Normal Can Reach

Returns \\\sup\|\gamma_1\|\\ over the skew normal family,
\$\$\dfrac{4-\pi}{2}\left(\dfrac{b}{\sqrt{1-b^2}}\right)^3 \approx
0.9952717,\$\$ with \\b = \sqrt{2/\pi}\\. The bound is attained only in
the limit \\\alpha \to \pm\infty\\, where the family degenerates to a
half-normal.

## Usage

``` r
sn_max_skew()
```

## Value

A single number.

## Details

A skewness beyond it belongs to no skew normal, so
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
bounds `gamma1` there and gives it a
[`linkfunctions7::bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html)
over the open interval. Without the bound the centered-to-direct map
would return a `NaN` several frames down, where the reason for it is no
longer visible.

The ceiling is approached slowly: measured, the skewness is 0.9556 at
\\\alpha = 10\\ and 0.99527 at \\\alpha = 10^4\\. Data skewer than this
needs
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

## See also

[`sn_b()`](https://statmodels7.github.io/distributions7/reference/sn_b.md)
for the constant it is built from,
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the parameter it bounds, and
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
for the family that goes further.

## Examples

``` r
distributions7:::sn_max_skew()
#> [1] 0.9952717

# The direct parametrization approaches it and does not pass it.
d1 <- skewnormal1_distrib()
vapply(c(10, 1e4, 1e8),
       function(a) skewness(d1, list(mu = 0, sigma = 1, alpha = a)), 0)
#> [1] 0.9555571 0.9952717 0.9952717

# It is where the centered parametrization's link is bounded.
skewnormal2_distrib()@params_bounds$gamma1
#> [1] -0.9952717  0.9952717
```
