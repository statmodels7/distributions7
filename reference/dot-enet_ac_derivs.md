# The Elastic Net in Its Two Rates

Returns the log-density's derivatives with respect to the location and
the two rates, \\(\mu, a, c)\\, up to the requested order. It is the
inner half of
[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md);
the outer half carries the result onto \\(\mu, \lambda, \alpha)\\.

## Usage

``` r
.enet_ac_derivs(y, p, order)
```

## Arguments

- y:

  A numeric vector of observations.

- p:

  The value of `.enet_parts()`, extended by
  [`.enet_g_higher()`](https://statmodels7.github.io/distributions7/reference/dot-enet_g_higher.md)
  when `order` is 3 or 4.

- order:

  A single integer, 1 to 4: the highest derivative order wanted.

## Value

A list of length `order`; element `k` is a named list of the order-`k`
derivatives, keyed as
[`deriv_names(c("mu", "a", "c"), k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them, each entry a numeric vector of the length of `y`.

## Why these coordinates

In them the log-density is \\-a\|z\| - cz^{2}/2 - \log 2 - L(x) +
\tfrac{1}{2}\log c\\ with \\z = y - \mu\\, \\x = ac^{-1/2}\\ and \\L =
\log M\\, so the data term is quadratic in \\z\\ and linear in each
rate. Every derivative of order three or more is therefore a derivative
of the **normalizing constant**, apart from
\\\partial^{3}\ell/\partial\mu^{2}\partial c = -1\\.

The normalizer is one pass of
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
over \\L(x(a,c))\\: the inner derivatives are \\-G, -G', -G'', -G'''\\
and the map has \\\partial^{2}x/\partial a^{2} = 0\\, so a key with two
or more \\a\\'s is an exact zero and the table leaves it out. The
remaining \\\tfrac12\log c\\ is pure in \\c\\ and is added component by
component.

## The kink

The location is not differentiable at \\z = 0\\, as in the Laplace, and
the sign function is what the first derivative carries there. Nothing at
order four touches \\\mu\\ at all.

## Notation

\\a = \lambda\alpha\\ is the Laplace rate, \\c = \lambda(1-\alpha)\\ the
Gaussian one, \\x = a/\sqrt c\\, \\L = \log M\\ with \\M\\ the Mills
ratio, and \\z = y - \mu\\.

## See also

[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md),
which calls this and carries the result onward;
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the partition sum; and
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
for the family.

## Examples

``` r
p <- distributions7:::.enet_g_higher(
  distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5)))
y <- c(-1.5, 0.4, 2.1)
out <- distributions7:::.enet_ac_derivs(y, p, 3L)
vapply(out, length, 0L)          # 3, 6 and 10 components
#> [1]  3  6 10

# The one third-order component that carries data.
out[[3]][["mu_mu_c"]]
#> [1] -1 -1 -1

# Everything else at order three is free of y, the data term being
# quadratic in z and linear in each rate.
vapply(out[[3]][setdiff(names(out[[3]]), "mu_mu_c")],
       function(v) diff(range(v)), 0)
#> mu_mu_mu  mu_mu_a   mu_a_a   mu_a_c   mu_c_c    a_a_a    a_a_c    a_c_c 
#>        0        0        0        0        0        0        0        0 
#>    c_c_c 
#>        0 
```
