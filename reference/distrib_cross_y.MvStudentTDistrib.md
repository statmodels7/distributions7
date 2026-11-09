# Multivariate Student t Mixed Response-Parameter Derivatives

Computes \\\partial^2\ell/\partial y\\\partial\theta_k\\, one \\n \times
p\\ matrix per parameter. The response gradient is \\-c\\w\\, so every
component carries the derivative of the weight beside the gaussian term
that it multiplies. With \\A_k = \partial\Sigma/\partial\eta_k\\,
\$\$\frac{\partial^2\ell}{\partial y\\\partial\mu_j} =
c\\\Sigma^{-1}e_j - \frac{\partial c}{\partial\mu_j}\\w, \qquad
\frac{\partial c}{\partial\mu_j} = \frac{2(\nu+p)\\w_j}{(\nu+q)^2},\$\$
\$\$\frac{\partial^2\ell}{\partial y\\\partial\eta_k} =
c\\\Sigma^{-1}A_k w - \frac{\partial c}{\partial\eta_k}\\w, \qquad
\frac{\partial c}{\partial\eta_k} = \frac{(\nu+p)\\w^\top A_k
w}{(\nu+q)^2},\$\$ \$\$\frac{\partial^2\ell}{\partial y\\\partial\nu} =
-\frac{(q-p)\\w}{(\nu+q)^2}.\$\$ For the gaussian the location block,
\\\Sigma^{-1}e_j\\, is the same for every row; here every component
varies with the row, and the location block carries the observation
through \\c\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters; each component is a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.
  The two differ here, because `nu` carries a log link by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\n \times p\\ numeric matrices, one per parameter, in
`distrib@params` order.

## Details

The log-density contains no distribution function, only `lgamma`, a
logarithm and a quadratic form, so every derivative in \\\nu\\ is closed
form. As \\\nu \to \infty\\ the weight tends to one and its derivatives
to zero, so every component that the two families share tends to the
gaussian's.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\\eta\\ the free vector of the matrix
parametrization, \\q\\ the squared Mahalanobis distance, \\c\\ the
weight, \\w = \Sigma^{-1}(y-\mu)\\ and \\e_j\\ the \\j\\th standard
basis vector.

## See also

[`distrib_grad_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvStudentTDistrib.md),
whose derivative in the parameters this is,
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md)
for the assembly,
[`distrib_cross_y.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvGaussianDistrib.md)
for the limiting family, and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

cy <- distrib_cross_y(d, y, theta)
names(cy)
#> [1] "mu1"          "mu2"          "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
#> [6] "nu"          

# The result is compared, parameter by parameter, with a difference of the
# response gradient.
h <- 1e-5
vapply(seq_along(d@params), function(k) {
  tp <- theta; tp[[k]] <- tp[[k]] + h
  tm <- theta; tm[[k]] <- tm[[k]] - h
  max(abs(cy[[k]] -
          (distrib_grad_y(d, y, tp) - distrib_grad_y(d, y, tm)) / (2 * h)))
}, numeric(1))
#> [1] 1.909717e-11 3.402589e-11 8.014878e-11 7.280865e-11 5.279199e-11
#> [6] 2.539313e-11

# As nu grows every shared component approaches the gaussian's.
t2 <- theta; t2$nu <- 1e7
g <- mvgaussian1_distrib(2)
cg <- distrib_cross_y(g, y, theta[1:5])
ct <- distrib_cross_y(d, y, t2)
max(vapply(1:5, function(k) max(abs(cg[[k]] - ct[[k]])), numeric(1)))
#> [1] 2.597118e-06
```
