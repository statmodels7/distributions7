# Does This Distribution Have a Real Quantile Method?

`TRUE` when the object gets its quantile function from a class-specific
method rather than from the numerical fallback.

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

Used to decide whether inverse transform sampling is cheap enough to
prefer over the ratio-of-uniforms sampler: inverting a *numerical*
quantile costs a call to `uniroot` per draw.

The test uses the documented S7 trick. S7 records the class a method was
registered on in the method's `signature` attribute, so an inherited
fallback is recognised by that class being `continuous_distrib` itself.
Note that [`identical()`](https://rdrr.io/r/base/identical.html) does
not work for this – S7 wraps the method object, so comparing against the
fallback fails even when it is the fallback.

## See also

[`rng_grou`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
