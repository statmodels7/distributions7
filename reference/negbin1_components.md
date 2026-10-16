# Derivative Components of NB1

Returns the components of
\\\partial^{a+b}\ell/\partial\mu^a\partial\theta^b\\ at any order from
one to four, from the sparse form the NB1 log-likelihood takes in the
size \\r = \mu/\theta\\.

`nb1_exact_cut()` returns the dispersion, 1, below which the
cancellation-free assembly is used in place of the recursion in the
size. The scores and second derivatives use the same threshold in the
compiled kernels.

`nb1_M_derivs()` returns \\M(\theta) = \log(1+\theta)/\theta\\ and its
derivatives to the given order. It is the one composite piece of the
cancellation-free form below, and it has a removable singularity at
\\\theta = 0\\: the recursion \\\theta M^{(b+1)} + (b+1)M^{(b)} = (-1)^b
b!/(1+\theta)^{b+1}\\ divides by \\\theta\\ and loses its digits there,
while the series \\M^{(b)} = (-1)^b \sum\_{n\ge 0} (-1)^n
\left\[\prod\_{i\le b}(n+i)\right\]\theta^n/(n+b+1)\\ converges only
below one. The function uses `nterm` terms of the series for \\\theta\\
below `cut` and the recursion otherwise.

`nb1_components_exact()` is the cancellation-free assembly, from the
form in \\\mu + i\theta\\ that `negbin1_psums_cpp()` documents.

## Usage

``` r
negbin1_components(y, theta, order)

nb1_exact_cut()

nb1_M_derivs(th, order, cut = 0.5, nterm = 80L)

nb1_components_exact(y, mu, th, order)
```

## Arguments

- y:

  A numeric vector of counts.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`, both strictly positive. Shorter
  components are recycled. Note that `theta` names both the list and its
  second component, the dispersion.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A named list of component vectors, one per distinct multi-index of the
given order and keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them: two at order 1, three at order 2, four at order 3 and five at
order 4. Each has the recycled length of the inputs.

## The sparse form, and the one composite piece

Writing \\r = \mu/\theta\\ the log-likelihood is \$\$\ell = G(r) + r
B(\theta) + C(\theta), \qquad G(r) = \log\Gamma(y+r) -
\log\Gamma(r),\$\$ with \\B = -\log(1+\theta)\\ and \\C = y\log\theta -
y\log(1+\theta)\\. The term \\rB(\theta)\\ is \\\mu B(\theta)/\theta\\,
linear in \\\mu\\, so it contributes to components carrying at most one
\\\mu\\; \\C\\ carries none. The only composite piece is
\\G(\mu/\theta)\\.

## The recursion in the powers of r

Its mixed derivatives take the form
\$\$\frac{\partial^{a+b}}{\partial\mu^a\partial\theta^b}G(\mu/\theta) =
\theta^{-(a+b)}\sum_j c_j\\ r^j\\ G^{(a+j)}(r),\$\$ and one further
\\\theta\\-derivative sends \$\$c_j r^j G^{(a+j)} \\\longrightarrow\\
-(a+b+j)\\c_j r^j G^{(a+j)} - c_j r^{j+1} G^{(a+j+1)},\$\$ because
\\\mathrm{d}r/\mathrm{d}\theta = -r/\theta\\ contributes to both the
power of \\r\\ and the order of \\G\\. The coefficients are integers,
and the recursion is run rather than solved, so every order is exact
with nothing transcribed beyond this one step.

## The cancellation the polygamma differences carry

Each \\G^{(m)}(r)\\ is a polygamma differenced at the shift \\y\\, which
is a count. As \\\theta \to 0\\ the family tends to the Poisson, \\r =
\mu/\theta\\ runs away, and the two terms of the difference agree to
leading order while the consumers above divide by \\\theta^{a+b}\\. The
differences therefore go through
[`psi_shift_diff()`](https://statmodels7.github.io/distributions7/reference/psi_shift_diff.md),
which forms them as an exact sum of reciprocals rather than as a
subtraction.

## The form without the size

The polygamma differences do not remove the cancellation among the
powers of \\r\\ in the recursion itself: at orders three and four those
terms are large compared with their sum, and as \\\theta\\ goes to zero
the recursion loses every significant digit.

The size need not appear at all, and the form without it has no such
cancellation. Since \\\log\Gamma(y+r) - \log\Gamma(r) =
\sum\_{i\<y}\log(r+i)\\ and \\r = \mu/\theta\\, the \\y\log\theta\\ each
such term carries cancels EXACTLY against \\C(\theta)\\, leaving
\$\$\ell = \sum\_{i\<y}\log(\mu + i\theta) - \log(y!) -
\frac{\mu}{\theta}\log(1+\theta) - y\log(1+\theta),\$\$ whose variable
does not run away and whose every derivative is elementary. Below
`nb1_exact_cut()` that is the route taken, through
`nb1_components_exact()`; above it the recursion is accurate and costs
\\O(1)\\ per observation, where the sum costs \\y\\ terms.

## See also

[`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
and
[`distrib_deriv4.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin1Distrib.md),
which call this;
[`psi_shift_diff()`](https://statmodels7.github.io/distributions7/reference/psi_shift_diff.md)
for the polygamma differences; and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the family.
