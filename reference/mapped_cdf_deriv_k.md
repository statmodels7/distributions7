# Higher CDF Derivatives of a Family Written as a Map of Another

The body every mapped family registers at orders 3 and 4, and at 2 where
its written-out route stops there. It asks whether the parent's cdf
derivatives are exact at **every** order up to the one wanted, takes
[`chain_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv_k.md)
if they are, and falls back to
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
on the new family's own cdf if they are not.

## Usage

``` r
mapped_cdf_deriv_k(
  distrib,
  parent,
  th_par,
  maps,
  q,
  theta,
  order,
  lower.tail,
  log,
  q_par = q
)
```

## Arguments

- distrib:

  The mapped family, whose cdf and parameter names are read.

- parent:

  The distribution it maps onto.

- th_par:

  The parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables.

- q:

  A numeric vector of quantiles, on the new family's own scale.

- theta:

  A named list of the new parameters.

- order:

  The derivative order, 2 to 4.

- lower.tail:

  Is the lower tail wanted? A single logical.

- log:

  Are derivatives of the log probability wanted? A single logical.

- q_par:

  The quantiles on the parent's scale. `q` by default, and `log(q)` for
  a lognormal.

## Value

A named list of numeric vectors of the requested order, on the requested
tail and scale.

## The gate

The chain is exact only if everything it carries is exact, and the
higher orders read the lower ones, so the test is over all of them. It
is what stops a family from adding an exact transformation to a
differenced parent and reporting the result as closed.

## A transformed response

A family may be the parent's law at a transformed point as well as at a
mapped parameter. A lognormal is a Gaussian at \\\log q\\, and since the
transformation carries no parameter, the derivatives in \\\theta\\ are
the parent's with the point substituted. `q_par` is where that
substitution is handed in.

## See also

[`chain_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv_k.md)
for the exact route;
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
for the fallback;
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
for the gate;
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md)
for orders 1 and 2.
