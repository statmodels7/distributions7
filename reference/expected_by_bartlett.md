# Expected Derivatives by the Bartlett Identity

Obtains the expected derivative of a given order from expectations of
*products of lower-order* derivatives, by summing over set partitions of
the index set.

## Usage

``` r
expected_by_bartlett(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations; only its length is used.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

## Value

A named list of expected derivative component vectors, each of length
`length(y)`.

## Details

The identity of order \\k\\ states that \$\$\sum\_{\pi}
\mathbb{E}\left\[\prod\_{B \in \pi} \ell_B\right\] = 0,\$\$ the sum
running over every partition \\\pi\\ of the index set. The single-block
partition is the target, so it equals minus the sum of all the others –
which is why
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
is the whole algorithm and why the top-order derivative is never needed.

At order 2 this reduces to \\\mathbb{E}\[\ell\_{ij}\] =
-\mathbb{E}\[\ell_i \ell_j\]\\, which needs only the score and remains
the information where the log-likelihood has a kink: there the pointwise
average of \\\partial^2 \ell\\ is not the information, while the score
variance still is.
[`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
calls this function at order 2 for such a family.

## See also

[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md),
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
