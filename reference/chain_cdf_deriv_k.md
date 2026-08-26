# The Chain Rule of Any Order on a Parent's CDF Derivatives

Carries a parent's derivatives of \\F\\ in its own parameters onto the
new ones, for a map given as keyed partial tables, at any order up to
four. The general form of
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md),
which writes orders 1 and 2 out; above them the same sum comes from
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md),
the enumeration the reparametrized parameter derivatives already run on,
so the package carries no second copy of the partition sum.

## Usage

``` r
chain_cdf_deriv_k(parent, q, th_par, maps, new_params, order)
```

## Arguments

- parent:

  The distribution being mapped.

- q:

  A numeric vector of quantiles, already on the parent's scale if the
  map transforms the response.

- th_par:

  The parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables. A missing key is an exact zero.

- new_params:

  A character vector naming the new parameters, which names and orders
  the result.

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(new_params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

The parent's tables are fetched on the natural scale, with
`lower.tail = TRUE` and `log = FALSE`, because the chain rule applies to
\\F\\ itself. The tail and the logarithm are put on afterwards by
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
so one Faa di Bruno pass serves both tails and both scales.

## Notation

\\F\\ is the parent's distribution function, \\\psi\\ the new
parameters, \\h\\ the map from them to the parent's, and \\\partial^I\\
a derivative with respect to a multi-index.

## See also

[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for orders 1 and 2;
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the partition sum;
[`mapped_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv_k.md),
the caller that gates it.
