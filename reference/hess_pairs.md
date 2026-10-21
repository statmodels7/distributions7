# Invert the Hessian Component Names

Maps each name produced by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
back to the pair of parameters with respect to which that component is
differentiated.

## Usage

``` r
hess_pairs(params)
```

## Arguments

- params:

  A character vector of parameter names.

## Value

A named list, parallel to `hess_names(params)`, in which each element is
a character pair.

## Details

Combining the parent's score with its Hessian in the wrappers requires
going from `"mu_sigma"` back to `c("mu", "sigma")`. Splitting the string
on `"_"` is the obvious way and it is wrong: a parameter whose own name
contains an underscore (`"log_scale"`) makes `"log_scale_log_scale"`
split into four pieces, and a function that takes the first and the last
piece returns `c("log", "scale")` without signaling an error. The map is
built from the parameter vector, so it does not depend on the characters
in the names.

[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
is the same idea for orders above two.

## See also

[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
