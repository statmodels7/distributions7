# Logistic Analytical Fourth-Order Derivatives

Closed-form observed fourth-order derivatives of the Logistic
log-density, in the notation of
[`distrib_deriv3.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md):
\$\$\dfrac{\partial^4 \ell}{\partial \mu^4} = \dfrac{g_4}{\sigma^4}\$\$
\$\$\dfrac{\partial^4 \ell}{\partial \mu^3 \partial \sigma} =
\dfrac{3g_3 + z g_4}{\sigma^4}\$\$ \$\$\dfrac{\partial^4 \ell}{\partial
\mu^2 \partial \sigma^2} = \dfrac{6g_2 + 6z g_3 + z^2 g_4}{\sigma^4}\$\$
\$\$\dfrac{\partial^4 \ell}{\partial \mu \partial \sigma^3} =
\dfrac{6g_1 + 18z g_2 + 9z^2 g_3 + z^3 g_4}{\sigma^4}\$\$
\$\$\dfrac{\partial^4 \ell}{\partial \sigma^4} = \dfrac{6 + 24z g_1 +
36z^2 g_2 + 12z^3 g_3 + z^4 g_4}{\sigma^4}\$\$

The expected fourth derivatives are not available in closed form, so
`expected = TRUE` is routed to
[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the approximated expected fourth
  derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
