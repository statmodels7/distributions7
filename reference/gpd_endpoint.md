# The Upper Endpoint of a Generalised Pareto

\\-\sigma/\xi\\ when \\\xi \< 0\\, and infinity otherwise.

## Usage

``` r
gpd_endpoint(sigma, xi)
```

## Arguments

- sigma:

  The scale, a positive numeric vector.

- xi:

  The shape, a numeric vector.

## Value

A numeric vector.

## Details

The endpoint depends on the parameters, which is the whole reason the
family needs care: for \\\xi \< 0\\ the support is bounded and moves
with \\\sigma\\ and \\\xi\\, so the licence to differentiate under the
integral sign is not automatic. See
[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
