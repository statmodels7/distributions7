# Derivatives of the Quadratic Form of a Multivariate Student t

Evaluates \\\partial^{B} q\\ for a multiset \\B\\ of mean and matrix
indices, with \\q = r'\Sigma^{-1}r\\. The form is quadratic in the mean,
so a block naming three or more mean coordinates is zero.

## Usage

``` r
mvt_q_deriv(b, r, pget, p)
```

## Arguments

- b:

  An integer vector of composite indices, mean coordinates first.

- r:

  The \\n \times p\\ matrix of residuals.

- pget:

  The accessor for \\\partial^{t}\Sigma^{-1}\\, as returned by
  [`mvg_ptensors`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md).

- p:

  The dimension.

## Value

A numeric vector of length `nrow(r)`, or a scalar recycled by the
caller.
