# Multivariate Student t Response Gradient

\\\partial \ell / \partial y = -c\\\Sigma^{-1}(y-\mu)\\, the gaussian
expression with the family's weight in front of it.

Closed form, from the family's own scale mixture rather than from a
sample.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  A numeric matrix of observations.

- theta:

  The parameters.

- scale:

  Either `"parameter"` or `"link"`.

- approx:

  Ignored; the answer is exact.

- nsim:

  Ignored; the answer is exact.

- ...:

  Ignored.

## Value

An \\n \times p\\ numeric matrix.

A named list of expected Hessian components.

## Details

Writing \\q = z^\top\Sigma^{-1}z\\ and \\w = (\nu+p)/(\nu+q)\\, the
variable \\v = q/(q+\nu)\\ is **exactly** \\\mathrm{Beta}(p/2, \nu/2)\\
and is independent of the direction \\z/\lVert z\rVert\\, which is
uniform on the sphere. Every expectation the information needs is
therefore a Beta moment or a polygamma, and the two averages separate:
the radial one gives \\wq = (\nu+p)v\\ and \\w^2q^2 = (\nu+p)^2v^2\\,
the angular one gives \\\operatorname{E}\[ee^\top\] = I/p\\ and
\\\operatorname{E}\[(e^\top Be)(e^\top Ce)\] =
\\\operatorname{tr}B\operatorname{tr}C + 2\operatorname{tr}(BC)\\/
\\p(p+2)\\\\.

With \\t_k = \operatorname{tr}(\Sigma^{-1}A_k)\\ and \\T\_{kl} =
\operatorname{tr}(\Sigma^{-1}A_k\Sigma^{-1}A_l)\\ the blocks come out as

\$\$I\_{\mu\mu} = \frac{\nu+p}{\nu+p+2}\\\Sigma^{-1}, \qquad I\_{kl} =
\frac{(\nu+p)T\_{kl} - t_kt_l}{2(\nu+p+2)}, \qquad I\_{k\nu} =
\frac{-t_k}{(\nu+p)(\nu+p+2)},\$\$

\$\$I\_{\nu\nu} = \frac{\psi'(\nu/2) - \psi'\\(\nu+p)/2\\}{4} -
\frac{p}{\nu(\nu+p)} + \frac{p}{2\nu(\nu+p+2)},\$\$

with the location orthogonal to everything else, every cross-expectation
with it being odd in \\z\\. Each reduces to the Gaussian's as \\\nu \to
\infty\\: \\I\_{\mu\mu} \to \Sigma^{-1}\\, \\I\_{kl} \to T\_{kl}/2\\,
and \\I\_{\nu\nu} \to 0\\, the degrees of freedom ceasing to be
identified in the limit.

## See also

[`mvstudent_t_distrib`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
