# Log-Tail Derivatives of the Location-Scale Families

Computes, in compiled code (`ls_log_tail_cpp()`), the derivatives of
order `o` of \\\log P\\ in \\(\mu, \theta_2)\\ for a family whose tail
probability is \\P = G(z)\\, \\z = (x - \mu)\\w(\theta_2)\\, with \\x =
q\\ or \\\log q\\. With \\g = \log G\\, every partial derivative is the
sum over the set partitions \\\pi\\ of the multi-index of
\\g^{(\|\pi\|)}(z)\prod\_{B \in \pi}\partial^B z\\; \\z\\ is linear in
\\\mu\\, so a block naming \\\mu\\ twice is zero.

## Usage

``` r
ls_log_tail_derivs(distrib, q, theta, o, spec, lower.tail)
```

## Arguments

- distrib:

  One of the families of `log_tail_ls`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with the location and the second parameter.

- o:

  The derivative order, 1 to 4.

- spec:

  The entry of `log_tail_ls` for the family.

- lower.tail:

  A single logical; `TRUE` selects the lower tail and `FALSE` the upper
  tail.

## Value

A named list, one component per multi-index of order `o`, in the order
of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order two.

## The bases

The lower tail of a symmetric base is its upper tail at \\-z\\. The
derivatives \\g^{(k)}\\, \\k \le 4\\, are written out:

- Laplace: \\g = \log\tfrac12 - z\\ for \\z \> 0\\ and \\\log(1 - u)\\,
  \\u = e^{z}/2\\, otherwise, whose derivatives are rational in \\u\\.

- logistic: polynomials in \\p = F\\ and \\s = 1 - F\\, where each is
  taken from its own tail by
  [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html).

- Gumbel, upper tail: with \\w = e^{-z}\\, \\g' = -\beta(w)\\,
  \\\beta(w) = w/(e^{w} - 1)\\, and each further order applies
  \\-w\\d/dw\\; below \\w = 1/2\\ by the Bernoulli series \\\sum_n n^k
  B_n w^n/n!\\, above it in closed form. Lower tail: \\\log F =
  -e^{-z}\\.

- Gaussian: below \\z = 9\\ from the hazard \\m = \phi/(1 - \Phi)\\ and
  \\m' = m(m - z)\\; from \\z = 9\\ by the asymptotic expansion \\g =
  -z^2/2 - \log z - \tfrac12\log 2\pi + \log\sum\_{n \le 40} (-1)^n
  (2n-1)!!\\z^{-2n}\\, whose truncation error there is below
  \\e^{-40}\\. For \\z \> 0\\ the term \\-z^2/2\\ is differentiated
  exactly as \\-(x - \mu)^2 w^2/2\\ and only the rest goes through the
  partition sum, which keeps the cancellation between the two out of the
  result.

Outside the support of a lognormal (\\q \le 0\\) the derivatives of
\\\log S\\ are zero and those of \\\log F\\ are `NaN`.

## Examples

``` r
d <- gaussian2_distrib()
# This is log S at 400 standard deviations, where the tail underflows.
distributions7:::ls_log_tail_derivs(d, 400, list(mu = 0, sigma2 = 1), 2L,
  list(base = "gaussian", w = "invsqrt"), lower.tail = FALSE)
#> $mu_mu
#> [1] -0.9999938
#> 
#> $sigma2_sigma2
#> [1] -160000.5
#> 
#> $mu_sigma2
#> [1] -400
#> 
```
