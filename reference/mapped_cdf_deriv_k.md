# Third and Fourth Log-CDF Derivatives of a Mapped Family

The chain rule on the parent's when the parent's are exact at every
order up to the one asked for, and the stencil otherwise.

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

  The mapped distribution.

- parent:

  The distribution being mapped.

- th_par:

  The parent's parameters at the new ones.

- maps:

  The map's keyed partial tables.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters.

- order:

  The derivative order, 3 or 4.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

- q_par:

  The point at which to evaluate the parent, when the two families are
  related by a monotone transformation of the response as well as by a
  map of the parameters. A lognormal is a gaussian at \\\log q\\, and
  since the transformation carries no parameter the derivatives in
  \\\theta\\ are the parent's with the point substituted.

## Value

A named list of derivative component vectors.

## Details

The gate is the one orders one and two use. A chain rule carrying a
differenced quantity would report a closed form and deliver the parent's
noise, and the truncation wrapper reads that distinction to choose its
own route.

## See also

[`mapped_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md)
