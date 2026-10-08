# The Derivatives of the Zero Wrappers' Expected Information

Compute the first derivatives of the expected information of
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
families in closed form, from the parent's expected information and its
derivatives and from the parent's mass, score and Hessian at zero, on
the parameter scale.

## Usage

``` r
zero_wrapper_dexpected(distrib, y, theta, threads = 1L)

zi_dexpected_entry(i, j, c, p, zi, f0, s, H, Ep, Dp)

za_disc_dexpected_entry(i, j, c, p, za, f0, s, H, Ep, Dp, n)

za_cont_dexpected_entry(i, j, c, p, za, Ep, Dp, n)
```

## Arguments

- distrib:

  A family from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  or
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  The response, read for its length.

- theta:

  The parameters, the parent's followed by the mixing probability.

- threads:

  The thread count passed to the parent's methods.

- i, j, c:

  The indices of one entry: the pair and the parameter differentiated,
  the parent's parameters first and the probability last.

- p:

  The number of the parent's parameters.

- zi, za:

  The wrapper's probability, one value per observation.

- f0:

  The parent's mass at zero.

- s, H:

  Functions of parameter indices returning the parent's score and
  Hessian at zero.

- Ep, Dp:

  Functions of parameter indices returning the parent's expected
  information and its derivatives.

- n:

  The number of observations.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
each component of length `length(y)`, on the parameter scale.

## Details

With \\f = f(0)\\, \\s\\ and \\H\\ the parent's score and Hessian at
zero, \\E\\ and \\D\\ the parent's expected information and its
derivatives, \\\zeta\\ the zero-inflation probability and \\L = \zeta +
(1-\zeta) f\\, the zero-inflated expected information is
\$\$E^{ZI}\_{ij} = (1-\zeta) E\_{ij} + \zeta(1-\zeta) f s_i s_j / L,\$\$
\$\$E^{ZI}\_{i\zeta} = -f s_i / L, \qquad E^{ZI}\_{\zeta\zeta} =
-(1-f)^2/L - (1-f)/(1-\zeta),\$\$ the observed Hessian at zero
cancelling from the first, and the derivatives follow with \\\partial_c
f = f s_c\\ and \\\partial_c s_i = H\_{ic}\\. For the zero-adjusted
discrete family, with \\q = 1 - f\\ and \\\pi\\ the probability of zero,
\$\$E^{ZA}\_{ij} = (1-\pi)\\\[E\_{ij}/q + f s_i s_j / q^2\], \qquad
E^{ZA}\_{\pi\pi} = -1/(\pi(1-\pi)),\$\$ the mixed block being zero; for
the continuous one \\E^{ZA}\_{ij} = (1-\pi) E\_{ij}\\.

Each entry is computed by `zi_dexpected_entry()`,
`za_disc_dexpected_entry()` or `za_cont_dexpected_entry()` with the
expressions of the scalar registry (`pt_wrappers.h`), so that the two
routes agree to the last bit.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
