# Derivative Arrays of a Dirichlet's Shape Vector

Collects \\\partial^{S}\alpha\\ for every multiset \\S\\ over the
Dirichlet's composite index set, the simplex's free values followed by
the concentration. The shape vector \\\alpha = \phi\\\mu(\eta)\\ is
BILINEAR in the concentration and the mean, so the whole array follows
from
[`simplex_map_tensors()`](https://statmodels7.github.io/distributions7/reference/simplex_map_tensors.md)
by three cases: \\\phi\\\partial^{S}\mu\\ when \\S\\ names no \\\phi\\,
\\\partial^{S'}\mu\\ when it names one, and exactly zero when it names
two or more.

## Usage

``` r
dirichlet_map_tensors(s, eta, phi, order)
```

## Arguments

- s:

  A
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  parametrization of the mean.

- eta:

  The mean's free vector, a numeric vector of length `s@n_free`.

- phi:

  The concentration, a single positive number.

- order:

  The highest order required, a single whole number from 1 to 4.

## Value

A named list of numeric vectors over the coordinates, keyed by sorted
tuple over the composite index set with \\\phi\\ LAST, so at
`s@n_free = 2` the concentration is index 3 and the key `"3,3"` holds a
vector of zeros.

## Notation

\\\alpha\\ is the Dirichlet's shape vector, \\\phi\\ its concentration,
\\\mu\\ the mean on the simplex, \\\eta\\ the mean's free vector, and
\\S'\\ the multiset \\S\\ with its one \\\phi\\ index removed.

## See also

[`simplex_map_tensors()`](https://statmodels7.github.io/distributions7/reference/simplex_map_tensors.md)
for the arrays it lifts,
[`tuple_indices_upto()`](https://statmodels7.github.io/distributions7/reference/tuple_indices_upto.md)
for the enumeration over the composite index set, and
[`dirichlet_higher()`](https://statmodels7.github.io/distributions7/reference/dirichlet_higher.md)
for the consumer.

## Examples

``` r
s <- parameters7::simplex(3)
eta <- c(0.3, -0.2)
dt <- distributions7:::dirichlet_map_tensors(s, eta, 8, 2)
names(dt)
#> [1] "1"   "2"   "3"   "1,1" "2,2" "3,3" "1,2" "1,3" "2,3"

# The key naming phi alone is the mean itself: alpha is phi times mu.
dt[["3"]]
#> [1] 0.4260125 0.2583897 0.3155978
as.numeric(parameters7::param_value(s, eta))
#> [1] 0.4260125 0.2583897 0.3155978

# Naming phi twice gives exactly zero, alpha being linear in it.
dt[["3,3"]]
#> [1] 0 0 0

# And a key naming no phi is phi times the mean's own derivative.
all.equal(dt[["1"]], 8 * as.numeric(parameters7::param_d1(s, eta)[[1]]))
#> [1] TRUE
```
