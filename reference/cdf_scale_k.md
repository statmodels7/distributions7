# CDF Derivatives on the Requested Tail and Scale, at Any Order

The general form of
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md):
converts derivatives of \\F\\ into derivatives of whichever tail was
asked for, on the natural or the logarithmic scale, at any order from 1
to 4.

## Usage

``` r
cdf_scale_k(distrib, Fq, dF, order, lower.tail, log)
```

## Arguments

- distrib:

  An object inheriting from `distrib`. Only its `params` are read, to
  name and enumerate the components.

- Fq:

  The distribution function at the quantile, a numeric vector.

- dF:

  A list of length `order`. Element \\k\\ is the table of \\k\\-th
  derivatives of \\F\\, keyed as
  [`deriv_names(distrib@params, k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
  on the natural scale and the lower tail.

- order:

  The derivative order wanted, 1 to 4.

- lower.tail:

  Is the lower tail wanted? A single logical.

- log:

  Are derivatives of the log probability wanted? A single logical.
  `FALSE` reads only `dF[[order]]` and the lower tables are then unused.

## Value

A named list of numeric vectors of the requested order alone. The lower
orders are consumed on the way and do not appear in the result.

## The two conversions

Switching to the upper tail flips the sign at every order, \\S = 1 -
F\\. Switching to the log scale is the moment-to-cumulant relation
\$\$\partial^I \log P = \sum\_\pi (-1)^{\|\pi\|-1}(\|\pi\|-1)! \prod\_{B
\in \pi} \frac{\partial^B P}{P},\$\$ summed over the set partitions of
the multi-index by
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md).
At second order it is the familiar \\\partial^2 P/P - (\partial
P/P)^2\\, and at third and fourth it has 4 and 15 terms.

## Why every lower order is needed

A partition into \\k\\ blocks reads \\k\\ ratios, so the relation at
order 4 needs the tables of orders 1, 2 and 3 as well. That is why `dF`
is a list of tables, and why a caller that wants only the fourth order
still assembles the first three.

## Notation

\\F\\ is the distribution function, \\S = 1 - F\\ the survival function,
\\P\\ whichever was asked for, \\\pi\\ a set partition of the
multi-index and \\B\\ one of its blocks.

## See also

[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md),
the written-out version for orders 1 and 2;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
for the partition sum;
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md),
which builds `dF`.
