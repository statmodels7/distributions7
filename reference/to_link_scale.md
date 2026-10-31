# Convert Parameter-Scale Derivatives to the Link Scale

Applies the multivariate Faa di Bruno formula with a diagonal Jacobian,
turning derivatives with respect to \\\theta\\ into derivatives with
respect to the unconstrained \\\eta\\.

## Usage

``` r
to_link_scale(distrib, theta, nat, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters, on the natural scale.

- nat:

  A list with `nat[[m]]` the named list of parameter-scale derivatives
  of order `m`. Its first element may be a list of zeros, as it is for
  expected derivatives.

- order:

  The derivative order to assemble, 1 to 6.

## Value

A named list of link-scale derivative component vectors.

## Details

The mathematics is set out in
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).
The implementation is not a plain transcription of the formula, for two
deliberate reasons.

First order is special-cased. It is a diagonal rescaling, \\\partial
\ell/\partial \eta_i = (\partial \ell/\partial \theta_i) h_i'\\, and it
is the order evaluated most often (once per scoring iteration), so it
bypasses the general assembly.

The index structure of the nested sum over \\j_t = 1, \dots, m_t\\
depends only on the parameter names and the order. It is computed once
by
[`link_scale_layout()`](https://statmodels7.github.io/distributions7/reference/link_scale_layout.md)
and cached, and the loop of this function performs only the arithmetic
on the vectors.

## See also

[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md),
[`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md)
