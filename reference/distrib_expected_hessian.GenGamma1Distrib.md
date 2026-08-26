# Generalized Gamma Expected Information

Returns the expected second derivatives in closed form. Every
expectation the observed Hessian needs is a moment of \\u = (Y/a)^{p}\\,
which is Gamma with shape \\k = d/p\\ and unit rate: \$\$E\[u\] = k,
\qquad E\[u\log u\] = k\\\psi(k+1), \qquad E\[u(\log u)^{2}\] =
k\\\psi(k+1)^{2} + \psi'(k+1)\\.\$\$ So `approx` and `nsim` are ignored
and
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
answers `TRUE`. It is the same device the Weibull and the Gumbel use,
where the corresponding variable is standard exponential and the
expectations are derivatives of \\\Gamma\\ at 2.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector. Its values do not enter the result, which is an
  expectation; only its length does, through recycling.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`. All three must be strictly
  positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  Ignored: the expectation is closed form. Accepted so that the
  signature matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of six numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `a_a`, `d_d`, `p_p`, `a_d`, `a_p`, `d_p`. Every entry is negative
definite as a matrix and free of the data.

## Identification

The three parameters are weakly identified together. \\d\\ and \\p\\
enter the density largely through their ratio \\k = d/p\\, and the
information reflects that: measured at \\a = 2, d = 3, p = 1.5\\ its
eigenvalues are 4.72, 0.138 and 0.00321, a condition number of 1470, and
the flat direction is \\(0.79, -0.52, 0.33)\\ in \\(a, d, p)\\. A fit of
all three wants several hundred observations; holding one with
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
is often the better model.

## Notation

\\k = d/p\\, \\\psi\\ and \\\psi'\\ the digamma and trigamma functions,
and \\u = (Y/a)^p\\ the Gamma-distributed transformation the family is
built on.

## See also

[`distrib_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md)
for the observed curvature,
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
for the predicate, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
th <- list(a = 2, d = 3, p = 1.5)
e <- distrib_expected_hessian(d, 0, th)
names(e)
#> [1] "a_a" "a_d" "a_p" "d_d" "d_p" "p_p"

# The scale-shape entry is the observed one, -1/a, the term being bilinear.
c(expected = e$a_d, observed = -1 / 2)
#> expected observed 
#>     -0.5     -0.5 

# The strategy argument is ignored, the expectation being closed form.
identical(e, distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 5))
#> [1] TRUE
distributions7:::expected_hessian_exact(d)
#> [1] TRUE

# The information is ill conditioned: d and p act largely through d/p.
nm <- c("a", "d", "p")
M <- matrix(c(e$a_a, e$a_d, e$a_p, e$a_d, e$d_d, e$d_p,
              e$a_p, e$d_p, e$p_p), 3, 3, dimnames = list(nm, nm))
ev <- eigen(-M)
c(eigenvalues = ev$values, condition = max(ev$values) / min(ev$values))
#> eigenvalues1 eigenvalues2 eigenvalues3    condition 
#> 4.721041e+00 1.379646e-01 3.211472e-03 1.470055e+03 
round(ev$vectors[, 3], 3)   # the flat direction
#> [1]  0.791 -0.517  0.327
```
