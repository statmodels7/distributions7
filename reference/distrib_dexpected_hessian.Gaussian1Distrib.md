# Gaussian Derivatives of the Expected Information

Returns the first and second derivatives of the expected information in
the parameters, in closed form from a compiled kernel. Only the standard
deviation moves it: \$\$\partial\_\sigma \mathbb{E}\[\ell\_{\mu\mu}\] =
\dfrac{2}{\sigma^3},\quad \partial\_\sigma
\mathbb{E}\[\ell\_{\sigma\sigma}\] = \dfrac{4}{\sigma^3},\quad
\partial\_{\sigma\sigma} \mathbb{E}\[\ell\_{\mu\mu}\] =
-\dfrac{6}{\sigma^4},\quad \partial\_{\sigma\sigma}
\mathbb{E}\[\ell\_{\sigma\sigma}\] = -\dfrac{12}{\sigma^4}.\$\$ Every
other component is zero. On the link scale the result is carried across
by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list with `mu` and `sigma`.

- scale:

  `"parameter"` or `"link"`.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused.

- threads:

  A single positive integer, the kernel's thread count.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
