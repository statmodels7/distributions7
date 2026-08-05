# The Largest Skewness a Skew Normal Can Reach

The supremum of \\\|\gamma_1\|\\ over the family, attained only in the
limit \\\alpha \to \pm\infty\\: \\(4-\pi)/2 \cdot (b/\sqrt{1-b^2})^3\\,
about 0.9952717.

## Usage

``` r
sn_max_skew()
```

## Value

A single number.

## Details

A skewness beyond it belongs to no skew normal, which is why the
constructor bounds the parameter there rather than letting the map
return a `NaN` several frames down. It is also the reason the skew \\t\\
exists.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
