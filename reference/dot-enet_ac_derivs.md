# The Elastic Net in Its Two Rates

The log-density's derivatives with respect to \\(\mu, a, c)\\ – the
location and the two rates – up to the requested order.

## Usage

``` r
.enet_ac_derivs(y, p, order)
```

## Arguments

- y:

  A numeric vector of observations.

- p:

  The value of `.enet_parts`, extended by `.enet_g_higher`.

- order:

  The derivative order, 1 to 4.

## Value

A list of length `order`; element `k` is the table of order-`k`
derivatives, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(c("mu", "a", "c"), k)`.

## Details

In these coordinates the log-density is \\-a\|z\| - cz^{2}/2 - \log 2 -
L(x) + \tfrac{1}{2}\log c\\ with \\z = y - \mu\\, \\x = ac^{-1/2}\\ and
\\L = \log M\\, so the data term is quadratic in \\z\\ and linear in
each rate. Every derivative of order three or more is therefore a
derivative of the normalizing constant, apart from
\\\partial^{3}\ell/\partial\mu^{2}\partial c = -1\\, and the normalizer
is one pass of
[`chain_assemble`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
over \\L(x(a,c))\\: the inner derivatives are \\-G, -G', -G'', -G'''\\
and the map has \\\partial^{2}x/\partial a^{2} = 0\\, which leaves a
short table.

The location is not differentiable at \\z = 0\\, as in the Laplace, and
the sign function is what the first derivative carries there.
