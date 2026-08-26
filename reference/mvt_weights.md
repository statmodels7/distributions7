# The Weight a Multivariate t Gives Each Observation

Computes the centered response \\r_i = y_i - \mu\\, its image \\w_i =
\Sigma^{-1}r_i\\ under the inverse scale matrix, the squared Mahalanobis
distance \\q_i = r_i^\top w_i\\, and the weight \\c_i = (\nu + p)/(\nu +
q_i)\\. Every derivative of the family is written in those four, so they
are formed once per call.

## Usage

``` r
mvt_weights(y, pc)
```

## Arguments

- y:

  An \\n \times p\\ numeric matrix of observations, already coerced by
  [`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md).

- pc:

  The result of
  [`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md),
  from which `mu`, `sigma_inv`, `nu` and `p` are read.

## Value

A named list with `r` and `w`, each an \\n \times p\\ numeric matrix,
and `q` and `cw`, each a numeric vector of length \\n\\.

## Details

The weight is the whole of the family's resistance to outliers. At \\q =
0\\ it is \\(\nu+p)/\nu\\, and it decays like \\1/q\\, so an observation
far from the location contributes less to every derivative instead of
dragging the fit towards itself. Letting \\\nu \to \infty\\ sends it to
one and recovers the gaussian, where nothing is downweighted.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\p\\ the dimension, \\q\\ the squared Mahalanobis
distance and \\c\\ the weight.

## See also

[`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md)
for the argument and
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
for the first consumer.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
              sigma_log_L2 = 0, sigma_L2.1 = 0, nu = 6)
pc <- distributions7:::mvt_pieces(d, theta)

# At the identity scale matrix q is the squared distance from the origin.
y <- rbind(c(0, 0), c(1, 0), c(3, 0), c(10, 0))
z <- distributions7:::mvt_weights(y, pc)
rbind(q = z$q, weight = z$cw)
#>            [,1]     [,2]      [,3]        [,4]
#> q      0.000000 1.000000 9.0000000 100.0000000
#> weight 1.333333 1.142857 0.5333333   0.0754717

# The weight starts at (nu + p) / nu and decays like 1 / q.
c(at_zero = (6 + 2) / 6, computed = z$cw[1])
#>  at_zero computed 
#> 1.333333 1.333333 

# As nu grows it flattens to one, and the family stops downweighting.
vapply(c(6, 60, 6000), function(nu) {
  t2 <- theta; t2$nu <- nu
  distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, t2))$cw[4]
}, numeric(1))
#> [1] 0.0754717 0.3875000 0.9839344
```
