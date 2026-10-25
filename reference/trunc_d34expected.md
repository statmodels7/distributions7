# The Third and Fourth Derivatives of a Truncated Family's Expected Information

Computes the derivatives of order `k` (3 or 4) of the expected Hessian
of a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family on the parameter scale, as expectations under the truncated law
taken by
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Usage

``` r
trunc_d34expected(d, y, th, k, t = 1L)
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

  The order, `3L` or `4L`.

- t:

  Unused; it is the thread count of the registration.

## Value

A list of components keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md),
each of length `length(y)`.

## Details

With \\\ell = \ell_P - \log Z\\ the truncated log-likelihood and
\\E\_{ab} = -E\[\ell_a\ell_b\]\\, the component for a multiset \\S\\ is
minus the expectation of the polynomial of
[`d34e_plan()`](https://statmodels7.github.io/distributions7/reference/d34e_plan.md),
which reads the derivatives of \\\ell\\ to order \\k + 1\\. The
truncation points do not depend on the parameters, so there is no
boundary term. Those of \\\ell\\ are the parent's less those of \\\log
Z\\, computed once per parameter combination by
[`trunc_logz_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_dexpected.md)
and passed to the integrand. The derivatives of \\\log Z\\ to order four
are read from the parent's distribution function where
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
holds, and the parent's distribution function has no derivatives of
order five, so each ratio \\d^B Z/Z\\ of that order is the expectation,
under the truncated law, of the complete Bell polynomial of the parent's
log-density derivatives (a quadrature or a sum), which is accurate but
not a closed form. Order five enters the fourth derivative only. The
form \\E\[\ell\_{ab}\]\\ would need the sixth derivative of the parent.

## See also

[`trunc_dexpected()`](https://statmodels7.github.io/distributions7/reference/trunc_dexpected.md)
