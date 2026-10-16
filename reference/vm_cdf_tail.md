# The von Mises Tails Below the Series' Resolution

`vm_cdf_tail()` returns \\P(Y \le q)\\ or \\P(Y \> q)\\, or their
logarithms, from the series of
[`vm_cdf()`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md)
where the requested tail is above \\10^{-6}\\, and from `vm_log_tail()`
below. The series has an absolute error of about \\10^{-16}\\, which
grows slowly with \\\kappa\\, so below that threshold its relative error
would exceed about \\10^{-9}\\. `vm_log_tail()` integrates
\\e^{\kappa\cos(y - \mu)}\\ over the tail, which excludes the mode, so
the integrand decreases from each end of the interval to the antimode.
The interval is split at the antimode when it contains it, and each
piece is integrated from its end on the scale \\1/\max\\\kappa\|\sin(e -
\mu)\|, \sqrt\kappa\\\\ at that end \\e\\, relative to the value there,
over panels of growing width until the remainder is below \\10^{-17}\\
of the sum. The difference of the cosines is written as a product of
sines. The pieces are added on the log scale, and the normalizing
constant comes from
[`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html),
so that the factors do not underflow at any concentration.

## Usage

``` r
vm_log_tail(q, mu, kappa, lower)

vm_log_piece(e, dir, L, mu, kappa)

vm_cdf_tail(q, mu, kappa, lower.tail, log.p)
```

## Arguments

- q:

  A numeric vector of angles in \\\[-\pi, \pi)\\.

- mu, kappa:

  Numeric vectors, the mean direction and the concentration.

- lower, lower.tail, log.p:

  Logical of length 1.

- e, dir, L:

  For `vm_log_piece()`, which integrates one piece: the end \\e\\ at
  which the integrand is largest, the direction (1 or -1) in which it
  decreases, and the length of the piece.

## Value

A numeric vector for `vm_cdf_tail()`; for `vm_log_tail()` and
`vm_log_piece()`, a number, the logarithm of the tail probability and of
the integral of \\e^{\kappa\cos(y - \mu)}\\ over the piece.

## Examples

``` r
distributions7:::vm_cdf_tail(-1, 0, 100, TRUE, FALSE)
#> [1] 5.101678e-22
```
