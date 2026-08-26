# One Univariate Chain Rule Over a Multivariate Map

Evaluates \\\partial^{S} f(u(v))\\ for a scalar function \\f\\ of one
coordinate, by the partition form of Faa di Bruno \$\$\partial^{S} f(u)
= \sum\_{\pi} f^{(\lvert\pi\rvert)}(u) \prod\_{B \in \pi} \partial^{B}
u,\$\$ the sum running over the set partitions of the multiset \\S\\.
Everything is vectorized over the coordinates, so one call serves every
\\j\\ at once.

## Usage

``` r
chain_univariate(tuple, fd, ud)
```

## Arguments

- tuple:

  An integer vector of free-value indices, with repeats allowed. Its
  length is the derivative order.

- fd:

  A list whose \\m\\-th element is \\f^{(m)}\\ evaluated at every
  coordinate, a numeric vector each. It must run to `length(tuple)`.

- ud:

  A named list of the map's derivative arrays, keyed by sorted tuple as
  [`simplex_map_tensors()`](https://statmodels7.github.io/distributions7/reference/simplex_map_tensors.md)
  returns them. Every block of every partition of `tuple` must have a
  key in it.

## Value

A numeric vector over the coordinates, as long as the elements of `fd`.

## Details

This is why the two simplex-valued families are cheap at orders three
and four. Their log-densities are sums of terms each depending on ONE
coordinate, \\\log\Gamma(\alpha_j)\\ for the Dirichlet and \\\log p_j\\
for the multinomial, so the multivariate Faa di Bruno collapses to one
univariate partition sum per coordinate and no mixed derivative array is
ever formed.

## Notation

\\S\\ is a multiset of free-value indices, \\\pi\\ a set partition of
it, \\B\\ a block of that partition, and \\u\\ the map being composed
with.

## See also

[`index_partitions()`](https://statmodels7.github.io/distributions7/reference/index_partitions.md)
for the enumeration it sums over,
[`simplex_map_tensors()`](https://statmodels7.github.io/distributions7/reference/simplex_map_tensors.md)
for the argument, and
[`dirichlet_higher()`](https://statmodels7.github.io/distributions7/reference/dirichlet_higher.md)
and
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
for the consumers.

## Examples

``` r
s <- parameters7::simplex(3)
eta <- c(0.3, -0.2)
ud <- distributions7:::simplex_map_tensors(s, eta, 2)

# Compose f(u) = log(u) with the map: f'(u) = 1/u, f''(u) = -1/u^2.
u <- as.numeric(parameters7::param_value(s, eta))
fd <- list(1 / u, -1 / u^2)

# The mixed second derivative of log(mu_j) in the two free values.
got <- distributions7:::chain_univariate(c(1L, 2L), fd, ud)
got
#> [1] 0.1100772 0.1100772 0.1100772

# Against a difference of the first-order chain rule, taken directly.
h <- 1e-5
d1 <- function(e) {
  uu <- as.numeric(parameters7::param_value(s, e))
  as.numeric(parameters7::param_d1(s, e)[[1]]) / uu
}
max(abs(got - (d1(eta + c(0, h)) - d1(eta - c(0, h))) / (2 * h)))
#> [1] 5.591458e-12
```
