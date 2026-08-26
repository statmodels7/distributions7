# Student t Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{2\nu}{\sigma^2(\nu+3)},\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \nu^2}\right\] =
\dfrac{1}{4}\left\[\psi_1\\\left(\dfrac{\nu+1}{2}\right) -
\psi_1\\\left(\dfrac{\nu}{2}\right)\right\] +
\dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma \\ \partial
\nu}\right\] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}.\$\$

The location is **orthogonal** to both other parameters: the two mixed
entries containing \\\mu\\ are exactly zero, the law being symmetric
about \\\mu\\ while the other two components are even in \\r\\. So
\\\hat\mu\\ is asymptotically independent of \\\hat\sigma\\ and
\\\hat\nu\\, which the scale and the degrees of freedom are not of each
other.

Because a closed form exists, `approx` and `nsim` are ignored: every
strategy returns the same six numbers. The arithmetic runs in a compiled
kernel, so the result does not depend on the thread count.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations. Only its length is read, the
  expectation not depending on the data; the values are ignored.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, the expectation being exact. Accepted so that the
  signature matches the generic's, where it selects between
  `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of six numeric vectors, `mu_mu`, `sigma_sigma`, `nu_nu`,
`mu_sigma`, `mu_nu` and `sigma_nu`, each of length `length(y)` and each
constant along it. `mu_sigma` and `mu_nu` are exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale, \\\nu \> 0\\ the degrees of freedom and
\\\psi_1\\ the trigamma function.

## Large degrees of freedom

The \\\nu\nu\\ entry is a difference of two trigammas half a unit apart
plus a term that cancels it to leading order, so writing it directly
loses every digit at a large \\\nu\\: before the asymptotic branches
were added it changed sign at \\\nu \approx 3\times10^5\\ and read
\\-n/2\\ on the link scale. The shipped kernel switches to a series at
measured crossovers and stays correctly signed to
`.Machine$double.xmax`.

## See also

[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the quantity this is the expectation of,
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)
for the score, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
eh <- distrib_expected_hessian(d, y, th)
vapply(eh, function(z) z[1], numeric(1))
#>        mu_mu  sigma_sigma        nu_nu     mu_sigma        mu_nu     sigma_nu 
#> -0.520833333 -0.868055556 -0.003022589  0.000000000  0.000000000  0.034722222 

# The four non-zero closed forms, written out.
c(mu_mu = -6 / (1.2^2 * 8),
  sigma_sigma = -10 / (1.2^2 * 8),
  nu_nu = (trigamma(3) - trigamma(2.5)) / 4 + 10 / (2 * 5 * 6 * 8),
  sigma_nu = 2 / (1.2 * 6 * 8))
#>        mu_mu  sigma_sigma        nu_nu     sigma_nu 
#> -0.520833333 -0.868055556 -0.003022589  0.034722222 

# Averaging the observed Hessian over draws reaches the same six numbers,
# the two containing mu going to zero.
set.seed(1)
z <- distrib_rng(d, 4e5, th)
vapply(distrib_hessian(d, z, th), mean, numeric(1))
#>         mu_mu   sigma_sigma         nu_nu      mu_sigma         mu_nu 
#> -5.211601e-01 -8.664218e-01 -3.030732e-03  1.067700e-03 -4.438806e-05 
#>      sigma_nu 
#>  3.461666e-02 

# The strategy argument is inert, the expectation being exact.
identical(eh, distrib_expected_hessian(d, y, th, approx = "mc", nsim = 50))
#> [1] TRUE
```
