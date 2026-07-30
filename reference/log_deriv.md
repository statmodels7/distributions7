# Derivatives of a Logarithm From the Ratios Alone

Computes \\d^I \log L\\ as \\\sum\_\pi (-1)^{\|\pi\|-1}(\|\pi\|-1)!
\prod\_{B \in \pi} (d^B L / L)\\.

## Usage

``` r
log_deriv(idx, ratio)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition.

- ratio:

  A function returning \\d^B L / L\\ for a block.

## Value

A numeric vector.

## Details

The moment-to-cumulant relation. What makes it the right tool here is
that only the **ratios** \\d^B L / L\\ are needed, never \\L\\'s own
derivatives – and the ratios are exactly what each wrapper can supply
cheaply, a truncated expectation for truncation and an affine expression
for the zero wrappers.

## See also

[`bell_f_ratio`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md),
the companion identity.
