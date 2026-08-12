# Elastic-Net Expected Information

The expected Hessian in closed form, which the family already carries
every piece of.

## Arguments

- distrib:

  An `EnetDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- scale:

  Either `"parameter"` or `"link"`.

- expected:

  Ignored; the answer is exact.

- approx:

  Ignored; the answer is exact.

- nsim:

  Ignored; the answer is exact.

- ...:

  Ignored.

## Value

A named list of expected Hessian components.

## Details

In the two rates the density is an exponential family with sufficient
statistics \\-\|z\|\\ and \\-z^2/2\\, so \\\log Z\\ is its cumulant
generating function and the information in \\(a, c)\\ is exactly the
Hessian of \\\log Z\\,

\$\$I\_{aa} = \operatorname{Var}(\|z\|), \quad I\_{ac} =
\operatorname{Cov}\\\left(-\|z\|, -\tfrac{z^2}{2}\right), \quad I\_{cc}
= \tfrac{1}{4}\operatorname{Var}(z^2),\$\$

which `.enet_logz_derivs()` computes for the observed Hessian already.
The map to \\(\lambda, \alpha)\\ is bilinear, so the information
transforms by \\J'IJ\\ with no second-derivative term.

Two of the three rate entries of the observed Hessian carry no data at
all and are therefore their own expectations, which is why
`lambda_lambda` and `alpha_alpha` below repeat them exactly. The third
does not: `lambda_alpha` carries \\-\|z\| + z^2/2\\, whose expectation
\\\partial_a \log Z - \partial_c \log Z\\ cancels the constant of the
same value sitting beside it and leaves the \\\log Z\\ term alone.

The location is where the two differ, and the reason is the kink. The
observed second derivative in \\\mu\\ is \\-c\\, which misses the point
mass \\\mathrm{d}\\\mathrm{sgn}(z)/\mathrm{d}z = 2\delta(z)\\ the
density carries at its own location, exactly as the Laplace does. The
information is defined there as the variance of the score, and

\$\$I\_{\mu\mu} = \operatorname{E}\left\[(a\\\mathrm{sgn}(z) +
cz)^2\right\] = a^2 + 2ac\\\operatorname{E}\|z\| +
c^2\operatorname{E}\[z^2\] = a^2 - 2ac\\\partial_a \log Z - 2c^2
\partial_c \log Z,\$\$

since \\\operatorname{E}\|z\| = -\partial_a \log Z\\ and
\\\operatorname{E}\[z^2\] = -2\partial_c \log Z\\. It reduces to
\\\lambda^2\\ at \\\alpha = 1\\, which is the Laplace's \\1/\sigma^2\\,
and to \\c\\ at \\\alpha = 0\\, which is the Gaussian's. The location
stays orthogonal to both rates, every cross-expectation vanishing by
symmetry.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
