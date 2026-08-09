# CDF Derivatives on the Requested Tail and Scale, at Any Order

The general form of
[`cdf_tail_scale`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md):
converts derivatives of \\F\\ of every order up to the one wanted into
derivatives of whichever tail was asked for, on the natural or the log
scale.

## Usage

``` r
cdf_scale_k(distrib, Fq, dF, order, lower.tail, log)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- Fq:

  The cdf evaluated at the quantile.

- dF:

  A list of length `order`; element `k` is the table of \\k\\-th
  derivatives of \\F\\, keyed as
  [`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(distrib@params, k)`.

- order:

  The derivative order wanted, 1 to 4.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

## Value

A named list of derivative component vectors of the requested order.

## Details

Switching tail flips the sign, since \\S = 1 - F\\. Switching to the log
scale is the moment-to-cumulant relation \\d^I \log P = \sum\_\pi
(-1)^{\|\pi\|-1}(\|\pi\|-1)! \prod\_{B} (d^B P / P)\\, which at second
order is the familiar \\d^2 P/P - (dP/P)(dP/P)\\ and at third and fourth
is what
[`log_deriv`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
sums. Only the ratios are needed, which is why every order up to `order`
has to be supplied.

## See also

[`cdf_tail_scale`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md),
[`log_deriv`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
