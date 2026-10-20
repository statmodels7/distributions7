# Derivative Arrays of a Simplex-Valued Map

Collects \\\partial^{S}\mu_j\\ for every multiset \\S\\ of free indices
up to the requested order, keyed by the sorted tuple, as a list of
numeric vectors over the coordinates \\j\\. The arrays come from
parameters7, which supplies orders three and four for a `simplex`. The
function derives nothing itself and only adds the keying, which lets the
chain rule in
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

  The highest order required, a single whole number from 1 to 4. A
  larger value signals an error.

## Value

A named list of numeric vectors, in which each vector is as long as the
simplex has coordinates. The names are `"1"`, `"2"`, `"1,1"`, `"1,2"`
and so on, and the indices in each name are sorted and comma-separated.
The length of the list is the number of index multisets of width 1 to
`order`.

## Notation

\\\mu\\ is the point of the simplex that the parametrization produces,
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

# Every derivative of a point on the simplex sums to zero, because the
# value sums to one at every eta.
vapply(mt, sum, numeric(1))
#>             1             2           1,1           2,2           1,2 
#>  2.775558e-17  2.775558e-17  6.938894e-18  1.387779e-17 -2.081668e-17 
```
