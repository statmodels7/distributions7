# The Elastic Net's Derivatives of a Given Order

Assembles the log-density's derivatives in \\(\mu, a, c)\\ through
[`.enet_ac_derivs()`](https://statmodels7.github.io/distributions7/reference/dot-enet_ac_derivs.md)
and carries them onto \\(\mu, \lambda, \alpha)\\ through
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
and the bilinear table of
[`.enet_rate_maps()`](https://statmodels7.github.io/distributions7/reference/dot-enet_rate_maps.md).
It is what
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
and
[`distrib_deriv4.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.EnetDistrib.md)
return.

## Usage

``` r
.enet_chain(y, theta, order)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, read
  positionally by `.enet_parts()`.

- order:

  A single integer, 1 to 4: the derivative order wanted.

## Value

A named list of the distinct order-`order` components in \\(\mu,
\lambda, \alpha)\\, named as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them, each a numeric vector of the length of `y`.

## Details

The license for orders three and four is that the **same assembly** run
at orders one and two reproduces the hand-written score and observed
Hessian: measured, exactly 0 at order one and \\1.1\times10^{-16}\\ at
order two. Against an independent route, one product stencil on the
analytic log-density, the pure-\\\lambda\\ components agree to
\\2\times10^{-4}\\ at order three and \\2\times10^{-3}\\ at order four,
which is the stencil's own accuracy at those orders.

## See also

[`.enet_ac_derivs()`](https://statmodels7.github.io/distributions7/reference/dot-enet_ac_derivs.md)
for the inner half,
[`.enet_rate_maps()`](https://statmodels7.github.io/distributions7/reference/dot-enet_rate_maps.md)
for the map,
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the partition sum, and
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
for the method that calls this.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# At order one it reproduces the hand-written score exactly, which is what
# licenses the orders where there is nothing to compare against.
g <- distrib_gradient(d, y, th)
a1 <- distributions7:::.enet_chain(y, th, 1L)
max(abs(unlist(a1[names(g)]) - unlist(g)))
#> [1] 0

# And at order two, the hand-written Hessian.
h <- distrib_hessian(d, y, th)
a2 <- distributions7:::.enet_chain(y, th, 2L)
max(abs(unlist(a2[names(h)]) - unlist(h)))
#> [1] 1.110223e-16
```
