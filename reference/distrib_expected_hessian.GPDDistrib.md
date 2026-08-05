# Generalized Pareto Analytical Expected Hessian

The closed form of Smith (1985), valid for \\\xi \> -1/2\\:
\$\$\mathbb{E}\left\[\dfrac{\partial^2\ell}{\partial\sigma^2}\right\] =
\dfrac{-1}{(1+2\xi)\sigma^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2\ell}{\partial\sigma\partial\xi}\right\]
= \dfrac{-1}{(1+2\xi)\sigma(1+\xi)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2\ell}{\partial\xi^2}\right\] =
\dfrac{-2}{(1+2\xi)(1+\xi)}\$\$ At \\\xi \le -1/2\\ the information does
not exist and `NA` is returned rather than a number; see
[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Arguments

- distrib:

  A `GPDDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `sigma` and `xi`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second-derivative components.

## See also

[`gpd_distrib`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
