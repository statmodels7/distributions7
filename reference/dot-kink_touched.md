# The Parameters with Which a Kink Moves

Returns the set \\\\p : \partial v/\partial \theta_p \ne 0\\\\, which it
evaluates at probe values inside each parameter's own interval and at a
few points of the support rather than assuming it.

## Usage

``` r
.kink_touched(distrib, kd)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- kd:

  A
  [`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md).

## Value

A character vector of parameter names, possibly empty.

## Details

The probe is the midpoint of a bounded interval and one unit inside a
half-line, which is the rule that
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
already uses; a component is called non-zero if it is non-zero at any
probe, since the set consists of the parameters for which the derivative
does not vanish identically.

## Examples

``` r
d <- laplace_distrib()
distributions7:::.kink_touched(d, kink_decomposition(d))
#> [1] "mu"
```
