# Expected Derivatives by Numerical Integration

Integrates each observed derivative component directly against the
density through
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
component by component. It is deterministic and is normally the most
accurate route where the observed derivative is available in closed
form.

## Usage

``` r
expected_by_integrate(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations. Only its length is used, to recycle
  the result: an expectation does not depend on the sample.

- theta:

  A named list of parameters, aligned to `distrib@params`.

- order:

  The derivative order, a single integer from 2 to 4.

## Value

A named list of expected derivative component vectors, each of length
`length(y)`, keyed as
[`observed_deriv()`](https://statmodels7.github.io/distributions7/reference/observed_deriv.md)
keys that order.

## Details

What this estimates is \\\mathbb{E}\[\partial^k \ell\]\\ literally. For
a regular model that is the quantity wanted. For a non-regular one it is
not the information, and the difference is total rather than small: on a
Laplace with no closed-form expected method, at \\\sigma = 1\\ over 200
observations, this returns **exactly 0** for the location component
while
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
returns \\-200 = -n/\sigma^2\\, which is the information. The observed
\\\ell\_{\mu\mu}\\ is zero almost everywhere, so its integral against
the density is zero and the point mass at the kink is invisible to it.

Quadrature is also unreliable where the observed derivative is itself a
finite difference, since it then integrates numerical noise. That error
is caught and re-raised naming the component and pointing at the
alternatives, so it reaches the caller as a statement about the route
rather than as an opaque failure from the integrator.

## See also

[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the three strategies;
[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md),
the route that survives a kink;
[`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
where quadrature struggles;
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which does the integration.
