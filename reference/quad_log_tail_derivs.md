# Log-Tail Derivatives by Quadrature of the Log-Density Derivatives

Computes the derivatives of order `o` of \\\log P\\, where \\P\\ is the
tail beyond `q` on the side that `lower.tail` selects, from \$\$\log P =
\ell(q) + \log J, \qquad J = \int\_{\mathrm{tail}} e^{D(y)}\\dy, \qquad
D(y) = \ell(y) - \ell(q).\$\$ The derivatives of \\\log J\\ follow from
the ratios \\\partial^B J / J\\, which are weighted means of the
complete Bell polynomials in the partial derivatives of \\D\\, through
the moment-to-cumulant sum; both sums are the recursions over the
subsets of a multi-index that contain its first position, run in
compiled code (`log_tail_cumulants_cpp()`) on the pairs that
`log_tail_pairs()` enumerates. Every partial of \\D\\ is a difference of
the family's own log-density derivatives at \\y\\ and at \\q\\, which
stays small over the tail, so the sums carry no cancellation of the size
that the ratios of
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
carry, and the computation does not underflow where \\P\\ does.

## Usage

``` r
quad_log_tail_derivs(distrib, q, theta, o, lower.tail, all_orders = FALSE)

quad_log_tail_prob(distrib, q, theta, lower.tail)
```

## Arguments

- distrib:

  A continuous distribution object.

- q:

  A numeric vector of quantiles inside the support.

- theta:

  A named list of parameters; each has the length of `q`.

- o:

  The derivative order, 1 to 4.

- lower.tail:

  A single logical; `TRUE` selects the lower tail and `FALSE` the upper
  tail.

- all_orders:

  Return the orders 1 to `o`, as a list of such lists?

## Value

A named list, one component per multi-index of order `o`, in the order
of
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
with the attribute `"lost"`.

`quad_log_tail_prob()` returns \\\log P = \ell(q) + \log J\\, a numeric
vector of the length of `q`.

## Details

The integral is taken with the 15-point Kronrod rule on 16 panels in \\u
\in \[0, 128\]\\. With \\s = 1/\|\partial_y \ell(q)\|\\ the scale on
which \\\ell\\ falls away from \\q\\, the map is \\y = q \pm s\\(e^{u} -
1)\\ toward an infinite end of the support, which turns an exponential
tail into a double exponential and a polynomial one into an exponential
in \\u\\, and \\y = b + (q - b)\exp\\-c\\(e^{u} - 1)\\\\, \\c = \min(1,
s/\|q - b\|)\\, toward a finite end \\b\\, which reaches \\b\\ at \\u =
\infty\\ and turns a power of \\y - b\\ into a double exponential. The
route is meant for the far tail, where the density decreases away from
`q`;
[`log_tail_route()`](https://statmodels7.github.io/distributions7/reference/log_tail_route.md)
uses it below \\\log P = -5\\. Toward an end at zero, a row whose last
node with a representable abscissa still carries a relative weight above
\\10^{-12}\\, while the nodes beyond it have underflowed to zero, has
mass below the smallest double (a shape near zero); it is marked in the
attribute `"lost"` of the result. Toward an infinite end the nodes past
the largest double are dropped, which leaves out a relative mass of
order \\(q/x\_{\max})^{\nu}\\ for a polynomial tail of index \\\nu\\.

`quad_tail_parts()` returns the nodes at every quantile, the response
vector \\(q, y_1, \dots)\\, in which the parameters are repeated to
match, the log densities and the weights \\e^{D}\\ times the Jacobian
and the rule's weight; `log_tail_nodes()` the nodes and weights in
\\u\\, computed once; `log_tail_pairs()` the components of orders 1 to
`o` in canonical form and, for each, the pairs (subset containing the
first position, rest) as indices into the components.

## Examples

``` r
d <- gamma1_distrib()
# This is log S of a gamma at 3000, where S underflows.
distributions7:::quad_log_tail_derivs(d, 3000, list(mu = 1.7, phi = 2.3),
                                      1L, lower.tail = FALSE)
#> $mu
#> [1] 451.6635
#> 
#> $phi
#> [1] 332.1425
#> 
#> attr(,"lost")
#> [1] FALSE
```
