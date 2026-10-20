# Standard Deviations and Correlations of a Structured Matrix

Turns a symmetric positive definite matrix and the derivatives of its
entries into standard deviations and correlations and computes the
Jacobian explicitly. Writing \\\Sigma = DRD\\ with \\D\\ the diagonal of
standard deviations, \$\$\frac{\partial s_j}{\partial\theta_k} =
\frac{A_k\[j,j\]}{2 s_j}, \qquad \frac{\partial
\rho\_{jl}}{\partial\theta_k} = \frac{A_k\[j,l\]}{s_j s_l} -
\frac{\rho\_{jl}}{2}\left(\frac{A_k\[j,j\]}{\Sigma\_{jj}} +
\frac{A_k\[l,l\]}{\Sigma\_{ll}}\right),\$\$ with \\A_k =
\partial\Sigma/\partial\theta_k\\. The function computes no further
derivative: the arrays are those that the caller obtained from the
parametrization.

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

  A \\p \times p\\ symmetric positive definite numeric matrix.

- a:

  A list of \\p \times p\\ numeric matrices, one per parameter of the
  distribution and in `params` order, holding
  \\\partial\Sigma/\partial\theta_k\\. For a parameter that does not
  affect the matrix, the entry is `NULL` or a matrix of zeros; both are
  accepted.

- params:

  The distribution's parameter names, which become the column names of
  the Jacobian.

- sd_label:

  The prefix for the diagonal quantities. Defaults to `"sd"`; the
  Student t passes `"scale_sd"`, because its diagonal quantities are not
  standard deviations of the response.

- cor_label:

  The prefix for the off-diagonal quantities. Defaults to `"cor"`.

- sd_block, cor_block:

  The block labels by which the printed summary groups the quantities.

## Value

A named list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
documents. There are \\p\\ diagonal quantities on the log scale and
\\p(p-1)/2\\ off-diagonal ones on Fisher's \\z\\.

## Notation

\\\Sigma\\ is the matrix, \\D\\ the diagonal of its square roots, \\R\\
the correlation matrix, \\s_j\\ a standard deviation, \\\rho\_{jl}\\ a
correlation, \\\theta\\ the parameter vector and \\A_k =
\partial\Sigma/\partial\theta_k\\.

## See also

[`mv_sigma_derivs()`](https://statmodels7.github.io/distributions7/reference/mv_sigma_derivs.md)
for the arrays that it takes,
[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
and
[`mv_derived.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvStudentTDistrib.md)
for the two callers, and
[`mv_entry_index()`](https://statmodels7.github.io/distributions7/reference/mv_entry_index.md)
for the labeling.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- distributions7:::align_theta(
  d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
          sigma_L2.1 = 0.5))
S <- mv_sigma(d, theta)
a <- distributions7:::mv_sigma_derivs(d, theta, 2)

sc <- distributions7:::mv_sd_cor(S, a, d@params)
sc$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.0000000 1.1180340 0.4472136 
sc$transform
#>     sd_v1     sd_v2 cor_v1_v2 
#>     "log"     "log"   "atanh" 

# The values really are the square roots of the diagonal and the
# correlations off it.
c(sqrt(diag(S)), S[2, 1] / sqrt(S[1, 1] * S[2, 2]))
#>        v1        v2           
#> 1.0000000 1.1180340 0.4472136 

# The Jacobian is compared with a numerical one taken through mv_sigma().
g <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  M <- mv_sigma(d, t2)
  c(sqrt(diag(M)), M[2, 1] / sqrt(M[1, 1] * M[2, 2]))
}
max(abs(sc$jacobian - numDeriv::jacobian(g, unlist(theta))))
#> [1] 8.833823e-12
```
