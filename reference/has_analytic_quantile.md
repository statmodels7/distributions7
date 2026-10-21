# Does This Distribution Have a Real Quantile Method?

Returns `TRUE` when the quantile function of the object comes from a
class-specific method rather than from the numerical fallback.

## Usage

``` r
has_analytic_quantile(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A single logical.

## Details

The result selects inverse transform sampling over the ratio-of-uniforms
sampler when the quantile is analytic, because inverting a *numerical*
quantile costs several evaluations of the distribution function per
draw, and each of them is a quadrature.

The test uses the documented S7 trick. S7 records the class on which a
method was registered in the method's `signature` attribute, so an
inherited fallback is recognized when that class is `continuous_distrib`
itself. The test does not use
[`identical()`](https://rdrr.io/r/base/identical.html), because S7 wraps
the method object, so comparing against the fallback fails even when the
method is the fallback.

## See also

[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
