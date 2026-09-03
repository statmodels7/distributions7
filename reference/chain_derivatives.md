# The Partition Sum Itself

Carries a parent's derivatives of the log-density into new coordinates,
given the keyed partial tables of the map. It is separated from
[`reparam_chain()`](https://statmodels7.github.io/distributions7/reference/reparam_chain.md)
so that a family written in its own right, not obtained through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
can use the same machinery instead of a second copy of it: the
hand-written second parametrizations all do.

## Usage

``` r
chain_derivatives(
  parent,
  y,
  th_par,
  maps,
  new_params,
  order,
  expected = FALSE,
  approx = "opg",
  nsim = 10000
)
```

## Arguments

- parent:

  The distribution whose derivatives are being carried.

- y:

  The response, a numeric vector.

- th_par:

  The parent's parameters as plain numbers, from
  [`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md).

- maps:

  The keyed partial tables of the map, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them.

- new_params:

  A character vector naming the new parameters, which names and orders
  the result.

- order:

  The derivative order, 1 to 4.

- expected:

  Should the expected derivatives be carried? A single logical, `FALSE`
  by default.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(new_params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## The identity

\$\$\ell^{(I)}(\psi) = \sum\_{\pi} \sum\_{i_1 \dots i\_{\|\pi\|}}
\ell^{(i_1 \dots i\_{\|\pi\|})}(\theta) \prod\_{B \in \pi}
\frac{\partial^{\|B\|}\theta\_{i_B}}{\partial\psi_B},\$\$ the sum over
set partitions of the multi-index. Every order the parent has in closed
form therefore survives in closed form.

## Under expectation

Expectation is linear and the map deterministic, so the expected
derivatives obey the same formula, and the term carrying the parent's
score drops because the score has mean zero. That is why `D[[1]]` is
left `NULL` when `expected` is `TRUE`: a parent with an exact expected
information gives an exact one here.

## Notation

\\\ell\\ is the log-density, \\\theta\\ the parent's parameters,
\\\psi\\ the new ones, \\\pi\\ a set partition of the multi-index and
\\B\\ one of its blocks.

## See also

[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md),
the sum it delegates to;
[`reparam_chain()`](https://statmodels7.github.io/distributions7/reference/reparam_chain.md),
its caller for a
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
object;
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md),
the same construction on the distribution function.
