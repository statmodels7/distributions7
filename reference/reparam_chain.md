# The Chain Rule of Any Order Through a Reparametrization

Carries a derivative of the parent's log-density into the new
coordinates.

## Usage

``` r
reparam_chain(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response.

- theta:

  A named list of the new parameters.

- order:

  The derivative order, 1 to 4.

- expected:

  Logical; if `TRUE`, carries the expected derivatives.

## Value

A named list of component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them.

## Details

With \\\theta = h(\psi)\\, the derivative of order \\\|I\|\\ is
\$\$\ell^{(I)}(\psi) = \sum\_{\pi} \sum\_{i_1 \dots i\_{\|\pi\|}}
\ell^{(i_1 \dots i\_{\|\pi\|})}(\theta) \prod\_{B \in \pi}
\frac{\partial^{\|B\|}\theta\_{i_B}}{\partial \psi_B}\$\$ the outer sum
running over the set partitions of the **positions** of \\I\\ and the
inner one over the assignment of a parent parameter to each block. This
is Faa di Bruno with a dense Jacobian, and it is the same partition
enumeration the wrappers of
[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and the rest already use.

Blocks index positions rather than variables, which is what makes a
repeated index carry its multiplicity without a factorial correction.

Expectation is linear and \\h\\ is deterministic, so the expected
derivatives obey the same formula with the parent's expected derivatives
in place of the observed ones. The first-order term drops, the score
having mean zero.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
