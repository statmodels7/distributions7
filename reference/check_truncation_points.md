# Validate the Truncation Endpoints

Checks that the interval is well formed and that truncating this parent
there leaves an identifiable model.

## Usage

``` r
check_truncation_points(distrib, lower, upper, is_disc)
```

## Arguments

- distrib:

  The parent distribution.

- lower:

  The lower endpoint, or `NULL`.

- upper:

  The upper endpoint, or `NULL`.

## Value

Invisibly `NULL`; raises an error on a bad interval.

## Details

The case worth naming is truncating zero away from a zero wrapper. The
\\(1-\zeta)\\ factor then cancels between the numerator and \\Z\\, so
\\\zeta\\ leaves the likelihood entirely and its score is identically
zero – the same defect as stacking the two zero wrappers, arriving by a
different route. Truncating a zero wrapper anywhere else is fine, and
the atom is carried through.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
