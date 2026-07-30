# Put CDF Derivatives on the Requested Tail and Scale

Converts derivatives of \\F\\ into derivatives of whichever tail was
asked for, on the natural or the log scale.

## Usage

``` r
cdf_tail_scale(distrib, Fq, dF1, dF2 = NULL, lower.tail, log)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- Fq:

  The cdf evaluated at the quantile.

- dF1:

  A named list of first derivatives of \\F\\.

- dF2:

  An optional named list of second derivatives of \\F\\.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

## Value

A named list of derivative component vectors.

## Details

Every route to a cdf derivative in this file produces derivatives of
\\F\\ itself; the `lower.tail` and `log` arguments are handled once,
here, rather than in each of them. Switching tail flips the sign, since
\\S = 1 - F\\, and switching to the log scale divides by the
probability, which for a second derivative brings in the familiar \\d^2
\log P = d^2 P / P - (dP/P)(dP/P)\\.

## See also

[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
