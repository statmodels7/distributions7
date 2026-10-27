# Student t Derivatives in the Standard Deviation

Return the derivatives of the log-density in \\(\mu, \sigma, \nu)\\ of
orders one to five, the expected information and its expected third and
fourth derivatives, and the first two derivatives of the expected
information, and each of them is computed by its own compiled kernel.

## Arguments

- distrib:

  A `StudentT2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  `"parameter"` or `"link"`; the link scale is applied by the generic.

- expected:

  Logical; for orders three and four, whether the expected derivative is
  returned.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list with one numeric vector per component: 3, 6, 10, 15 and 21
components at orders one to five. At order two the components are keyed
`mu_mu`, `sigma_sigma`, `nu_nu`, `mu_sigma`, `mu_nu`, `sigma_nu`. The
expected information has the same six components, the expected third and
fourth derivatives have 10 and 15,
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
has 18 and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
has 36.

## Details

With \\z = (y - \mu)/\sigma\\, \\k = \nu - 2\\ and \\q = z^2/k\\, the
log-density is \$\$\ell = c(\nu) - \log\sigma - \tfrac12\log\pi -
\frac{k + 3}{2}\log(1 + q), \qquad c(\nu) =
\log\Gamma\\\left(\frac{\nu+1}{2}\right) -
\log\Gamma\\\left(\frac{\nu}{2}\right) - \tfrac12\log(\nu - 2).\$\$
Every component is a closed form derived offline and written out per
component, one kernel per order. The part in the data is a rational
function of \\z\\ and \\k\\, reduced symbolically and written in \\q\\,
\\1/k\\ and \\t = 1/(1 + q)\\; the score in \\\nu\\ carries the
logarithm through \\q/(1+q) - \log(1 + q)\\. The derivatives of
\\c(\nu)\\, and the expected derivatives in \\\nu\\, are evaluated from
the polygamma functions below \\\nu = 20\\ and from their asymptotic
series in \\1/\nu\\ above it, because their terms cancel to leading
order as \\\nu\\ grows. The expectations are closed forms: under the
model \\t\\ follows a beta distribution with parameters \\\nu/2\\ and
\\1/2\\.

## See also

[`distrib_pdf.StudentT2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentT2Distrib.md).

## Examples

``` r
d <- student_t2_distrib()
th <- list(mu = 1, sigma = 2, nu = 6)
y <- c(-2, 0.5, 4)
g <- distrib_gradient(d, y, th)
h <- 1e-6
up <- distrib_pdf(d, y, list(mu = 1, sigma = 2, nu = 6 + h), log = TRUE)
dn <- distrib_pdf(d, y, list(mu = 1, sigma = 2, nu = 6 - h), log = TRUE)
all.equal(g$nu, (up - dn) / (2 * h), tolerance = 1e-6)
#> [1] TRUE
```
