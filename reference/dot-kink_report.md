# Print the Table of check_kink

Prints one line per check, with the worst gap behind it.

## Usage

``` r
.kink_report(out, gaps, jumps, pred, kd, ord, distrib)
```

## Arguments

- out:

  The logical vector of results.

- gaps:

  Per-parameter gaps of the `dv` check.

- jumps:

  Measured jumps of the score, or `NULL`.

- pred:

  Predicted jumps, or `NULL`.

- kd:

  A
  [`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md).

- ord:

  The order the composition leaves.

- distrib:

  The family.

## Value

`NULL`, invisibly; called for the printing.

## Examples

``` r
check_kink(laplace_distrib(), verbose = TRUE)
#> check_kink: laplace,  phi = abs,  order 0
#>   dv against a difference of v      [PASSED]  worst 4.55e-12
#>   the jump of the score             [PASSED]  worst 0.00e+00
#>   the smooth parameters             [PASSED]
#>       mu         measured            2   declared            2
#>       sigma      measured            0   declared            0
```
