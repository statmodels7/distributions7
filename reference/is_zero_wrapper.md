# Does This Distribution Already Model a Probability of Zero?

`TRUE` for a distribution produced by
[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
or
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
in either of its two forms.

## Usage

``` r
is_zero_wrapper(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A single logical.

## See also

[`check_not_stacked`](https://statmodels7.github.io/distributions7/reference/check_not_stacked.md)
