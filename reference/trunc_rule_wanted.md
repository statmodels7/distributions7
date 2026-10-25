# Does a Truncated Family Take Its Expectations on the Compiled Rule?

Returns `TRUE` for a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
continuous family whose parent is a skew t or a Student t of class
`StudentT1Distrib`, the families whose tail is \\\|z\|^{-1-\nu}\\ and
for which the compiled rule carries anchor nodes.

## Usage

``` r
trunc_rule_wanted(distrib)
```

## Arguments

- distrib:

  A truncated family.

## Value

A single logical.
