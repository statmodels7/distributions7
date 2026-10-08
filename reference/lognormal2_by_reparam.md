# Lognormal Distribution in the Mean and Variance of Y, Obtained

The same family as
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md),
obtained through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
on
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
rather than written out.

## Usage

``` r
lognormal2_by_reparam(link_mean = log_link(), link_var = log_link())
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_var:

  Link function for the variance. Defaults to the log.

## Value

A reparametrized distribution object.

## Details

This exists as a reference for the tests and is not exported:
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
carries its own kernels, so the two are independent implementations of
one law.

## See also

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
