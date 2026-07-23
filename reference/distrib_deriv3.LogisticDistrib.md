# Logistic Analytical Third-Order Derivatives

Closed-form observed third-order derivatives of the Logistic
log-density.

Writing \\z = (y-\mu)/\sigma\\, \\t = 1/(1+e^{-z})\\ and \\u = 1-t\\,
the log-density is \\\ell = -\log\sigma + g(z)\\ with \\g(z) = -z -
2\log(1+e^{-z})\\, whose derivatives in \\z\\ are \$\$g_1 = 1-2t, \quad
g_2 = -2tu, \quad g_3 = -2tu(1-2t), \quad g_4 = -2tu(1-6tu).\$\$ Since
\\z\\ depends on both parameters, each component is a polynomial in
\\z\\ with the \\g_j\\ as coefficients: \$\$\dfrac{\partial^3
\ell}{\partial \mu^3} = -\dfrac{g_3}{\sigma^3}\$\$ \$\$\dfrac{\partial^3
\ell}{\partial \mu^2 \partial \sigma} = -\dfrac{2g_2 + z
g_3}{\sigma^3}\$\$ \$\$\dfrac{\partial^3 \ell}{\partial \mu \partial
\sigma^2} = -\dfrac{2g_1 + 4z g_2 + z^2 g_3}{\sigma^3}\$\$
\$\$\dfrac{\partial^3 \ell}{\partial \sigma^3} = -\dfrac{2 + 6z g_1 +
6z^2 g_2 + z^3 g_3}{\sigma^3}\$\$

The expected third derivatives are not available in closed form, so
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

  Logical; if `TRUE`, returns the approximated expected third
  derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
