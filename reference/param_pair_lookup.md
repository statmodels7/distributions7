# Where Each Pair of Free Values Sits in a Parametrization's Second Derivatives

Builds a lookup from a pair of free-value positions to the position of
the matching component of
[`parameters7::param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.html).
A Hessian method walks the unordered pairs of the distribution's own
parameters and has to find, for each matrix pair, the second derivative
array that belongs to it; the two enumerations are not the same, because
the distribution's parameters carry the mean components in front.

## Usage

``` r
param_pair_lookup(s)
```

## Arguments

- s:

  A parameters7 parametrization, from which `param_tuple_indices()`
  supplies the enumeration.

## Value

A named list of single integers, one per unordered pair of free values,
keyed `"k:l"` with \\k \le l\\. Its length is
`s@n_free * (s@n_free + 1) / 2`.

## Details

The keys are built from parameters7's own enumeration rather than by
taking a component name apart, for the reason that package documents: a
free value whose label contains the separator splits into the wrong
number of pieces, so a name is not a safe route back to a pair of
indices.

## See also

[`mv_hess_indices()`](https://statmodels7.github.io/distributions7/reference/mv_hess_indices.md)
for the other half of the same bookkeeping and
[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the consumer.

## Examples

``` r
s <- parameters7::log_cholesky(2)
lk <- distributions7:::param_pair_lookup(s)
unlist(lk)
#> 1:1 2:2 3:3 1:2 1:3 2:3 
#>   1   2   3   4   5   6 

# The lookup is a permutation of the positions of param_d2()'s components.
d2 <- parameters7::param_d2(s, c(0.1, -0.2, 0.4))
setequal(unlist(lk), seq_along(d2))
#> [1] TRUE

# The array it finds for the pair (1, 3) is the mixed second derivative in
# those two free values, against a difference of the first derivatives.
eta <- c(0.1, -0.2, 0.4)
h <- 1e-5
num <- (parameters7::param_d1(s, eta + c(0, 0, h))[[1]] -
        parameters7::param_d1(s, eta - c(0, 0, h))[[1]]) / (2 * h)
max(abs(d2[[lk[["1:3"]]]] - num))
#> [1] 2.655653e-13
```
