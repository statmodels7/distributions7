# Derivatives of the Covariance with Respect to Every Parameter

Returns \\\partial\Sigma/\partial\theta_k\\ for each parameter of a
multivariate distribution built on a covstructs7 structure, as a list
aligned with `distrib@params` and `NULL` where the covariance does not
depend on the parameter.

## Usage

``` r
mv_sigma_derivs(distrib, theta, n_before)
```

## Arguments

- distrib:

  A distribution carrying a `struct` property.

- theta:

  A named list of parameters, already aligned.

- n_before:

  How many parameters precede the structure's free values.

## Value

A list of matrices and `NULL`s, of length `distrib@n_params`.

## Details

The mean components and, for a Student \\t\\, the degrees of freedom
leave the matrix alone, so those entries are `NULL` and cost nothing.
When the structure parametrises the precision the chain rule of an
inverse applies, \\\partial\Sigma/\partial\eta_k = -\Sigma A_k \Sigma\\.
