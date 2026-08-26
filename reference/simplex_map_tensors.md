# Derivative Arrays of a Simplex-Valued Map

Collects \\\partial^{S}\mu_j\\ for every multiset \\S\\ of free indices
up to the requested order, keyed by the sorted tuple, as a list of
numeric vectors over the coordinates \\j\\. The arrays come from
parameters7, whose contract carries orders three and four for a
`simplex`, so nothing is derived here; what the function adds is the
keying, which lets the chain rule in
[`chain_univariate()`](https://statmodels7.github.io/distributions7/reference/chain_univariate.md)
look a block up by its index multiset.

## Usage

``` r
simplex_map_tensors(s, eta, order)
```

## Arguments

- s:

  A
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  parametrization.

- eta:

  The free vector, a numeric vector of length `s@n_free`.

- order:

  The highest order required, a single whole number from 1 to 4. Above 4
  the `switch` returns `NULL` and the call fails.

## Value

A named list of numeric vectors, each as long as the simplex has
coordinates, keyed `"1"`, `"2"`, `"1,1"`, `"1,2"` and so on, with the
indices sorted and comma-separated. Its length is the number of index
multisets of width 1 to `order`.

## Notation

\\\mu\\ is the point of the simplex the parametrization produces,
\\\eta\\ its free vector, and \\S\\ a multiset of free-value indices.

## See also

[`chain_univariate()`](https://statmodels7.github.io/distributions7/reference/chain_univariate.md),
which consumes these,
[`dirichlet_map_tensors()`](https://statmodels7.github.io/distributions7/reference/dirichlet_map_tensors.md)
for the version that carries a concentration as well, and
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
for the parametrization.

## Examples

``` r
s <- parameters7::simplex(3)
eta <- c(0.3, -0.2)
mt <- distributions7:::simplex_map_tensors(s, eta, 2)
names(mt)
#> [1] "1"   "2"   "1,1" "2,2" "1,2"

# The first-order entry is the parametrization's own first derivative.
mt[["1"]]
#> [1]  0.2445259 -0.1100772 -0.1344486
as.numeric(parameters7::param_d1(s, eta)[[1]])
#> [1]  0.2445259 -0.1100772 -0.1344486

# Every derivative of a point on the simplex sums to zero, the value
# summing to one at every eta.
vapply(mt, sum, numeric(1))
#>             1             2           1,1           2,2           1,2 
#>  2.775558e-17  2.775558e-17  6.938894e-18  1.387779e-17 -2.081668e-17 
```
