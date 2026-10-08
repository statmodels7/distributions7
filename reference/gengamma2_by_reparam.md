# Generalized Gamma Distribution in the Mean, Obtained

The same family as
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md),
obtained through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
on
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
rather than written out.

## Usage

``` r
gengamma2_by_reparam(
  link_mean = log_link(),
  link_d = log_link(),
  link_p = log_link()
)
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_d:

  Link function for the shape. Defaults to the log.

- link_p:

  Link function for the power. Defaults to the log.

## Value

A reparametrized distribution object.

## Details

This exists as a reference for the tests and is not exported:
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md)
carries its own kernels, so the two are independent implementations of
one law.

## See also

[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md)
