# Logistic Analytical Observed Hessian

Computes the analytical observed Hessian (second derivatives) of the
Logistic log-density with respect to the parameters \\\mu\\ and
\\\sigma\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2}
\text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right)\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2}
\left\[ 1 - \dfrac{(y-\mu)^2}{2\sigma^2}
\text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) -
\dfrac{2(y-\mu)}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right)
\right\]\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} =
-\dfrac{1}{\sigma^2} \left\[ \tanh\left(\dfrac{y-\mu}{2\sigma}\right) +
\dfrac{y-\mu}{2\sigma} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right)
\right\]\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of second derivatives.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
