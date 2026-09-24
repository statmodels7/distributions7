# Derivatives of the Expected Information, Elementary Families

The first and second derivatives of the expected information in the
parameters, from compiled kernels, for seventeen families whose expected
information is an elementary function of the parameters. Each component
is an ordinary derivative of the family's written-out
\\\mathbb{E}\[\ell\_{ab}\]\\; on the link scale the result is carried
across by
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## Arguments

- distrib:

  A distribution object of one of the classes above.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

- threads:

  A single positive integer, the kernel's thread count.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

The expected informations differentiated, on the parameter scale:

- Bernoulli and binomial of size \\n\\: \\-n/(\mu(1-\mu))\\, with \\n =
  1\\ for the Bernoulli;

- exponential: \\-1/\mu^2\\; geometric: \\-1/(\mu(1+\mu))\\;

- chi-squared: \\-\psi'(\mu/2)/4\\;

- Cauchy, logistic and Gumbel: every component a constant over
  \\\sigma^2\\, so only \\\sigma\\ moves them;

- gaussian by its variance, and lognormal in \\(\mu, \sigma^2)\\:
  \\-1/\sigma^2\\ and \\-1/(2\sigma^4)\\; gaussian by its precision:
  \\-\tau\\ and \\-1/(2\tau^2)\\;

- inverse gaussian: \\-1/(\phi\mu^3)\\ and \\-1/(2\phi^2)\\ by the
  dispersion, \\-\lambda/\mu^3\\ and \\-1/(2\lambda^2)\\ by the shape;

- gamma by its mean and variance: with \\a = \mu^2/\sigma^2\\ and \\g(a)
  = a\\\psi'(a) - 1\\, \\\mathbb{E}\[\ell\_{\mu\mu}\] =
  -(1+4g)/\sigma^2\\, \\\mathbb{E}\[\ell\_{\mu\sigma^2}\] = 2\mu
  g/\sigma^4\\ and \\\mathbb{E}\[\ell\_{\sigma^2\sigma^2}\] = -a
  g/\sigma^4\\; \\g\\ and its derivatives are formed from the remainders
  of the polygamma functions' asymptotic series, so nothing cancels at a
  large shape, where \\g \sim 1/(2a)\\;

- beta by its two shapes: \\\psi'(\alpha+\beta) - \psi'(\alpha)\\,
  \\\psi'(\alpha+\beta)\\ and \\\psi'(\alpha+\beta) - \psi'(\beta)\\;

- Weibull by its scale and shape: \\-\sigma^2/\mu^2\\,
  \\(1-\gamma)/\mu\\ and \\-C/\sigma^2\\, with \\\gamma\\ Euler's
  constant and \\C = (1-\gamma)^2 + \pi^2/6\\;

- generalized Pareto by its scale and shape: with \\d = 1 + 2\xi\\,
  \\-1/(d\sigma^2)\\, \\-1/(d\sigma(1+\xi))\\ and \\-2/(d(1+\xi))\\, and
  `NA` for \\\xi \le -1/2\\, where the information does not exist.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
