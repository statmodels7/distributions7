# Multivariate Student t Score

Computes the first derivatives of the log-density in closed form. With
\\w = \Sigma^{-1}(y-\mu)\\, \\q = (y-\mu)^\top w\\ and \\c =
(\nu+p)/(\nu+q)\\, \$\$\partial\_\mu \ell = c\\w, \qquad
\partial\_{\eta_k}\ell =
-\tfrac{1}{2}\partial\_{\eta_k}\log\lvert\Sigma\rvert + \tfrac{c}{2}\\
w^\top A_k w,\$\$ \$\$\partial\_\nu \ell = \tfrac{1}{2}\left\[
\psi\\\left(\tfrac{\nu+p}{2}\right) -
\psi\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{\nu} -
\log\\\left(1+\tfrac{q}{\nu}\right) +
\tfrac{(\nu+p)q}{\nu(\nu+q)}\right\].\$\$ The gaussian's score is the
limit \\c \to 1\\. Every observation enters the location and matrix
components through that one weight, and its decay with \\q\\ is the
whole of the family's resistance to a distant row.

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

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here: under the default log link the
  `nu` component is multiplied by \\\nu\\, and the rest are unchanged.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector of length \\n\\ per parameter, in
`distrib@params` order: \\p\\ location components, one per free value of
the matrix parametrization, then `nu`.

## Details

The \\\nu\\ component is NOT evaluated in the form printed above. Both
of its bracketed pairs cancel to leading order as \\\nu\\ grows, so it
is assembled as \\A_p(\nu) + D(q/\nu) + (p/\nu)u/(1+u)\\ instead, with
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
and
[`mvt_D()`](https://statmodels7.github.io/distributions7/reference/mvt_D.md)
carrying out each cancellation exactly. The direct form loses every
digit at large \\\nu\\ and changes sign at \\\nu = 10^9\\; see
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
for the measurements.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\\eta\\ the free vector of the matrix
parametrization, \\A_k = \partial\Sigma/\partial\eta_k\\, \\q\\ the
squared Mahalanobis distance, \\c\\ the weight, \\\psi\\ the digamma
function and \\\ell\\ the log-density of one observation.

## See also

[`distrib_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md)
for the observed curvature,
[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md)
for the information,
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the weight,
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
for the cancellation, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 25, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>          mu1          mu2 sigma_log_L1 sigma_log_L2   sigma_L2.1           nu 
#>    5.8944688    0.5141700   -2.8238795   -9.5932528    2.2949745    0.2040964 

# Against a numerical derivative of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 8.141114e-09

# The location component is the weighted whitened residual, so a distant
# row contributes less than a near one of the same direction.
z <- distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, theta))
all.equal(cbind(g$mu1, g$mu2), z$cw * z$w, check.attributes = FALSE)
#> [1] TRUE

# The nu component decays as nu^-2, computed without cancellation.
vapply(c(1e3, 1e6, 1e10), function(nu) {
  t2 <- theta; t2$nu <- nu
  distrib_gradient(d, y[1, ], t2)$nu
}, numeric(1))
#> [1] 3.393052e-07 3.394278e-13 3.394279e-21

# And the link scale differs, this family having one link that is not the
# identity.
c(parameter = sum(g$nu),
  link = sum(distrib_gradient(d, y, theta, scale = "link")$nu))
#> parameter      link 
#> 0.2040964 1.2245781 
```
