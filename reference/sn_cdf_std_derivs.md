# The Skew Normal's Distribution Function in Standard Coordinates

Returns \\\partial^i_z \partial^j\_\alpha G\\ for \\G(z, \alpha) =
\Phi(z) - 2T(z, \alpha)\\ over the pairs with \\1 \le i + j \le\\
`order`, with \\T\\ Owen's T. This is the whole of the skew normal's cdf
derivative surface in standard coordinates; the location and the scale
are chained on afterwards by
[`sn_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_deriv_k.md).

## Usage

``` r
sn_cdf_std_derivs(z, al, order)
```

## Arguments

- z:

  The standardized quantile, a numeric vector.

- al:

  The shape, a numeric vector of any sign recyclable against `z`.

- order:

  The highest total order wanted, 1 to 4.

## Value

A nested list indexed `[[i + 1]][[j + 1]]`, holding
\\\partial^i_z\partial^j\_\alpha G\\ as a numeric vector. Entries with
\\i + j\\ above `order` are absent.

## Why the integral never has to be differentiated

Owen's T is defined by an integral, and both of its partial derivatives
are elementary: \$\$\frac{\partial T}{\partial h} =
-\varphi(h)\left\\\Phi(ah) - \tfrac12\right\\, \qquad \frac{\partial
T}{\partial a} = \frac{\varphi(h)\\\varphi(ah)}{1+a^2}.\$\$ The integral
is therefore differentiated away at the first order and never has to be
differentiated again. A quantity defined by an integral is not thereby a
quantity whose derivatives need one.

## The first derivatives, and the check on them

\\G_z = 2\varphi(z)\Phi(u)\\ and \\G\_\alpha = -2\varphi(z)\varphi(u)R\\
with \\u = \alpha z\\ and \\R = (1+\alpha^2)^{-1}\\. The first is the
skew normal's own density, which confirms the identity above; the second
is checked against the gradient the family reports.

Everything above the first order is one Leibniz rule over those factors,
with \\\Phi(u)\\ and \\\varphi(u)\\ differentiated by Faa di Bruno over
\\u\\. That map is bilinear, so it has only three non-zero partials and
the expansion stays short.

## Notation

\\z\\ is the standardized quantile, \\\alpha\\ the shape, \\u = \alpha
z\\, \\\Phi\\ and \\\varphi\\ the standard normal distribution and
density, and \\T\\ Owen's T.

## See also

[`sn_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_deriv_k.md),
which chains this onto the location and the scale;
[`recip_1p_sq()`](https://statmodels7.github.io/distributions7/reference/recip_1p_sq.md)
for the shape factor;
[`numericals7::owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.html)
for the function itself.
