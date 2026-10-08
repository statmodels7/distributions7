# The Derivatives of a Truncated Family's Expected Information

Computes the first (`k = 1`) or second (`k = 2`) derivatives of the
expected information of a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family on the parameter scale, as expectations under the truncated law
taken by
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
the sums over the support or the quadratures that give the expected
information itself.

## Usage

``` r
trunc_dexpected(d, y, th, k, t = 1L)

trunc_logz_derivs(d, theta, order)
```

## Arguments

- d:

  A truncated family, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  The response, read for its length.

- th:

  The parameters.

- k:

  The order, `1L` or `2L`.

- t:

  Unused; the thread count of the registration.

- theta:

  The parameters.

- order:

  The highest order, 1 to 4.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
(`k = 1`) or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
(`k = 2`), each component of length `length(y)`.

## Details

With \\\ell\\ the truncated log-likelihood, \\\ell - \log Z\\, and \\G =
\ell_i \ell_j\\, the second Bartlett identity gives \\E\[\ell\_{ij}\] =
-E\[G\]\\, and differentiating under the expectation, \$\$\partial_c
E\[G\] = E\[G_c + G\ell_c\],\$\$ \$\$\partial\_{cd} E\[G\] =
E\[G\_{cd} + G_c\ell_d + G_d\ell_c + G(\ell_c\ell_d + \ell\_{cd})\],\$\$
with \\G_c = \ell\_{ic}\ell_j + \ell_i\ell\_{jc}\\ and \\G\_{cd} =
\ell\_{icd}\ell_j + \ell\_{ic}\ell\_{jd} + \ell\_{id}\ell\_{jc} +
\ell_i\ell\_{jcd}\\. The truncation points do not depend on the
parameters, so no boundary term arises. The truncated derivatives are
the parent's less those of \\\log Z\\, which do not depend on \\y\\:
they are computed once per parameter combination by
`trunc_logz_derivs()` and passed to the integrand, which reads only the
parent's derivatives at the nodes. Registered through
[`register_dexpected()`](https://statmodels7.github.io/distributions7/reference/register_dexpected.md).

`trunc_logz_derivs()` returns the derivatives of \\\log Z\\ to order
`order`, keyed by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
at each order (the pairs as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
spells them), from the ratios \\d^B Z / Z\\ by
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md):
a ratio is read from the parent's distribution function where its
derivatives of that order exist, and is otherwise the expectation of the
complete Bell polynomial under the truncated law.

## See also

[`distrib_expected_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedContinuousDistrib.md)
