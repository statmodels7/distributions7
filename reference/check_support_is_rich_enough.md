# Reject a Model With More Parameters Than the Support Can Distinguish

Enforces the counting rule that makes a zero wrapper identifiable.

## Usage

``` r
check_support_is_rich_enough(distrib, fun)
```

## Arguments

- distrib:

  The parent distribution being wrapped.

- fun:

  The calling constructor's name, used in the message.

## Value

Invisibly `NULL`; raises an error when the support is too small.

## Details

A discrete distribution on \\k\\ points has \\k-1\\ free probabilities,
and either wrapper spends `n_params + 1` of them, so \\k \ge\\
`n_params + 2` is necessary. The bound is the same for inflation and for
adjustment.

What it rules out is exactly the Bernoulli, and
`binomial_distrib(size = 1)` with it. Zero-inflating a Bernoulli gives
two parameters for the one free cell of \\\\0, 1\\\\; zero-adjusting it
leaves the truncated part concentrated on \\\\1\\\\ with no free
parameter at all, so \\\mu\\ disappears from the likelihood and the pmf
is literally the same for \\\mu = 0.2\\ and \\\mu = 0.9\\.

## See also

[`n_support_points`](https://statmodels7.github.io/distributions7/reference/n_support_points.md),
[`check_not_stacked`](https://statmodels7.github.io/distributions7/reference/check_not_stacked.md)
