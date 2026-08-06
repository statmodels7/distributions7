# Beta Mixed Derivatives in Mean and Precision

Closed form. The shapes are \\a = \mu\phi\\ and \\b = (1-\mu)\phi\\, and
\\\ell^{(y)} = (a-1)/y - (b-1)/(1-y)\\, so the components are \\\phi/y +
\phi/(1-y)\\ and \\\mu/y - (1-\mu)/(1-y)\\.

## Arguments

- distrib:

  A `Beta1Distrib` object.

- y:

  A numeric vector in the unit interval.

- theta:

  A list containing `mu` and `phi`.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
