# Derivative Components of the Generalized Pareto

Returns the components of
\\\partial^{a+b}\ell/\partial\sigma^a\partial\xi^b\\ at any order from
one to four, by splitting the log-density into a part that is analytic
at \\\xi = 0\\ and a part that carries a removable singularity there.

## Usage

``` r
gpd_components(y, theta, order, cut = 0.2, threads = 1L)
```

## Arguments

- y:

  A numeric vector of observations on the support.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`. `sigma` must be strictly
  positive; `xi` may be of either sign, including zero, which is the
  exponential limit. Shorter components are recycled.

- order:

  The derivative order, an integer from 1 to 4.

- cut:

  The value of \\\lvert\xi z\rvert\\ below which the series is used,
  defaulting to `0.2`. It is exposed so that a test can force either
  route in the region where both are accurate and compare them; a caller
  has no reason to change it.

- threads:

  A single positive integer, how many threads the polynomial kernel may
  use. Defaults to `1L`.

## Value

A named list of component vectors, one per distinct multi-index of the
given order and keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them: two at order 1, three at order 2, four at order 3 and five at
order 4. Each has the recycled length of the inputs.

## The split

With \\z = y/\sigma\\ and \\t = 1 + \xi z\\, \$\$\ell = -\log\sigma -
L - W, \qquad L = \log t, \qquad W = L/\xi.\$\$ \\t\\ is affine in
\\\xi\\, so every partial of it carrying two or more \\\xi\\ vanishes
and the written-out template of
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
covers \\L\\ outright.

## Why W has two routes

\\W = L/\xi\\ has a removable singularity at \\\xi = 0\\, where it tends
to \\z\\. Differentiating it by Leibniz, \$\$\partial^{a,b} W =
\sum\_{j} \binom{b}{j}\\ \partial^{a,j}L \cdot
\partial^{b-j}(1/\xi),\$\$ produces terms of size \\\xi^{-(b+1)}\\ that
cancel against each other, so near zero the form loses every digit.
Below a threshold it is replaced by the series the cancellation leaves,
\$\$W = \sum\_{k \ge 0} \frac{(-1)^k \xi^k z^{k+1}}{k+1},\$\$
differentiated term by term.

The threshold is on \\\lvert\xi z\rvert\\ and **not** on \\\xi\\. At
order \\b\\ the Leibniz terms are of size \\z\xi^{-b}\\ against an
answer of size \\z^{b+1}\\, so the relative cancellation is \\(\xi
z)^{-b}\\ and the accuracy is \\\varepsilon(\xi z)^{-b}\\, which at
order four reaches 0.3 by \\\xi z = 1.7\times10^{-4}\\. That prediction
was measured, and matching it is what fixed the switch. The series needs
\\\lvert\xi z\rvert \< 1\\ to converge and is cut at 0.2, where forty
terms leave \\10^{-28}\\ while the Leibniz form still carries
\\10^{-13}\\, so the two overlap.

## What the series branch evaluates

The two elementwise powers the sum appears to need are algebraically
one: with \\u = \xi z\\, \\\xi^{k-b}z^{k+1} = u^{k-b}z^{b+1}\\, so
\\z^{b+1}/\sigma^a\\ leaves the sum and what remains is a **polynomial
in** \\u\\ whose coefficients are scalar in \\k\\. `gpd_poly_cpp()`
evaluates it by Horner from the highest power down, which sums a
decaying series smallest first. Each branch also runs on its own
elements only, so a sample straddling the cut does not pay for both in
full.

## See also

[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
and
[`distrib_deriv4.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GPDDistrib.md),
which call this;
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
for the two-variable composition template; and
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for the family.
