# Standard Deviations and Correlations of a Structured Matrix

Turns a symmetric positive definite matrix and the derivatives of its
entries into the standard deviations and correlations it decomposes
into, together with their Jacobian.

## Usage

``` r
mv_sd_cor(
  sigma,
  a,
  params,
  sd_label = "sd",
  cor_label = "cor",
  sd_block = "Standard deviations",
  cor_block = "Correlations"
)
```

## Arguments

- sigma:

  The matrix.

- a:

  A list of derivative matrices, one per parameter, in the order of the
  distribution's parameters. Entries may be `NULL` for parameters the
  matrix does not depend on.

- params:

  The parameter names, used to label the Jacobian's columns.

- sd_label:

  The label the diagonal quantities are named with.

- cor_label:

  The label the off-diagonal quantities are named with.

- sd_block, cor_block:

  The headings the two groups print under.

## Value

A list as described in
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).

## Details

Writing \\\Sigma = D R D\\ with \\D\\ the diagonal of standard
deviations, the two readings are \\s_j = \sqrt{\Sigma\_{jj}}\\ and
\\\rho\_{jk} = \Sigma\_{jk}/(s_j s_k)\\. Differentiating,
\$\$\frac{\partial s_j}{\partial \eta_k} = \frac{A_k\[j,j\]}{2 s_j},
\qquad \frac{\partial \rho\_{jk}}{\partial \eta_l} =
\frac{A_l\[j,k\]}{s_j s_k} - \frac{\rho\_{jk}}{2}
\left(\frac{A_l\[j,j\]}{\Sigma\_{jj}} +
\frac{A_l\[k,k\]}{\Sigma\_{kk}}\right),\$\$ with \\A_k =
\partial\Sigma/\partial\eta_k\\.
