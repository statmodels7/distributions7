# Poisson Derivatives of the Expected Information

The first and second derivatives of \\\mathbb{E}\[\ell\_{\mu\mu}\] =
-1/\mu\\, \\1/\mu^2\\ and \\-2/\mu^3\\, from a compiled kernel. On the
link scale the result is carried across by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A `PoissonDistrib` object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list with `mu`.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

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
