# The Order Left by a Composition

The largest `k` for which a composition has `k` continuous derivatives
at the origin: 0 for \\\lvert v \rvert\\ and \\(v)\_+\\, whose first
derivative jumps; 1 for \\(v)\_+^2\\; 2 for \\(v)\_+^3\\; and -1 for
\\1\\v\>0\\\\, which is not continuous at all.

## Usage

``` r
kink_order(phi)
```

## Arguments

- phi:

  A single string naming the composition; see
  [`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md).

## Value

A single number, `-1`, `0`, `1` or `2`.

## Details

The order -1 is worse in kind than the others. There a smoothing does
not repair a missing derivative but a discontinuity of the log-density
itself, and the unsmoothed problem has no well-defined maximum, so it is
reported apart rather than treated with the rest.

## See also

[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md),
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md).

## Examples

``` r
kink_order("abs")
#> [1] 0
kink_order("hinge3")
#> [1] 2
```
