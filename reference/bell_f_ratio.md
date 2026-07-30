# Complete Bell Polynomial in the Parent's Log-Derivatives

Computes \\d^I f / f\\ from the derivatives of \\\log f\\, as the sum
over set partitions \\\sum\_\pi \prod\_{B \in \pi} \ell^{(B)}\\.

## Usage

``` r
bell_f_ratio(idx, ell)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition.

- ell:

  A function returning the parent's log-derivative for a block.

## Value

A numeric vector.

## Details

This is the Bartlett lemma read backwards: instead of using the identity
to eliminate a derivative, it uses it to build one. Every wrapper needs
it, because each of their log-likelihoods is the parent's log-density
plus, or instead of, \\\log L\\ for some \\\theta\\-dependent \\L\\.

## See also

[`log_deriv`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the companion identity.
