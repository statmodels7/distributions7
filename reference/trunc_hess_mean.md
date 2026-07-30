# Truncated Mean of the Parent's Hessian

\\\mathbb{E}\_T\[H\_{ij}\]\\ for every Hessian component, by quadrature.

## Usage

``` r
trunc_hess_mean(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- theta:

  A named list of parameters.

## Value

A named list, one component per Hessian entry.
