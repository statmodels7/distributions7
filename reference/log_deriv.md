# Derivatives of a Logarithm From the Ratios Alone

Computes \\d^I \log L\\ from the ratios \\d^B L / L\\, as \$\$d^I \log L
= \sum\_{\pi} (-1)^{\|\pi\|-1}(\|\pi\|-1)! \prod\_{B \in \pi} \frac{d^B
L}{L},\$\$ the moment-to-cumulant relation. Only the **ratios** are
needed, never \\L\\'s own derivatives and never \\L\\ itself, so a
wrapper can differentiate a normalizing constant it can only evaluate up
to scale.

## Usage

``` r
log_deriv(idx, ratio)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition, naming the
  multi-index \\I\\. Its length is the order.

- ratio:

  A function of one block, returning \\d^B L / L\\ for that block as a
  numeric vector. Called once per block of every partition.

## Value

A numeric vector, the length of whatever `ratio` returns: \\d^I \log L\\
at each observation.

## Details

It is the inverse of
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md):
that sum carries derivatives of a logarithm to derivatives of the
function, this one carries them back. At order 1 it returns \\d_i L /
L\\ and at order 2 \\d\_{ij}L/L - (d_i L/L)(d_j L/L)\\, the familiar
relation for the second derivative of a logarithm.

## See also

[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md),
the companion identity;
[`index_partitions()`](https://statmodels7.github.io/distributions7/reference/index_partitions.md),
which supplies the partitions;
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
the wrappers that consume it.
