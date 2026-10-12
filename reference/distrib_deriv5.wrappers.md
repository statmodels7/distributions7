# The Wrappers' Fifth Derivatives

Computes the fifth derivatives of the log-likelihood of
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
and
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
families. The zero wrappers, the truncated and the folded families use
the partition sums that give their third and fourth derivatives, one
order up, and carry them to the link scale with
[`deriv5_scale()`](https://statmodels7.github.io/distributions7/reference/deriv5_scale.md);
a transformed and a fixed family read the parent's fifth derivative
directly, on the scale requested, since their parameters and links are
the parent's.

## Arguments

- builder:

  One of the order-generic builders of `R/wrapper_derivatives.R` or
  [`fold_deriv_k()`](https://statmodels7.github.io/distributions7/reference/fold_deriv_k.md).

- distrib:

  A family from one of the wrappers above.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of fifth-derivative components keyed by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
at order 5, each of length `length(y)`.

## Details

The zero wrappers, the truncated and the folded families assemble the
fifth derivative of \\\log L\\ from the block ratios \\d^B L / L\\ by
the moment-to-cumulant relation over the 52 partitions of five indices,
each ratio a complete Bell polynomial in the parent's derivatives to
order five. A transformed family's fifth derivatives are the parent's at
the preimage, and a fixed family's are the parent's at the full
parameter vector, subset to the free parameters. The result is analytic
where the parent's fifth derivative is.

The partition sums are used only where the parent's fifth derivative is
analytic, as
[`has_exact_deriv5()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv5.md)
reports. Where it is the stencil of
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
the products of the Bell polynomials and, for a truncated family, the
quadrature of its ratios amplify the stencil's error, and the wrapper's
fifth derivative is instead one stencil on its own analytic fourth.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
