# Inverse Gaussian Distribution in the Mean and Shape, Obtained

The same family as
[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md),
obtained through
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
rather than written out.

## Usage

``` r
invgauss2_by_reparam()
```

## Value

A reparametrized distribution object.

## Details

This exists as a check rather than as a second way of doing the same
thing.
[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
carries its own kernels, so the two are independent implementations of
one object and their agreement needs no tolerance to be chosen. It is
not exported for that reason.

## See also

[`invgauss2_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
