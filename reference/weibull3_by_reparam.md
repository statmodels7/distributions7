# Weibull Distribution in the Mean, Obtained

The same family as
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md),
obtained through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
on
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
rather than written out.

## Usage

``` r
weibull3_by_reparam(link_mean = log_link(), link_sigma = log_link())
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_sigma:

  Link function for the shape. Defaults to the log.

## Value

A reparametrized distribution object.

## Details

This exists as a reference for the tests and is not exported:
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
carries its own kernels, so the two are independent implementations of
one law.

## See also

[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
