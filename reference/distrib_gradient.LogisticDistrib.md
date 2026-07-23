# Logistic Analytical Gradient

Computes the analytical gradient (first derivatives) of the Logistic
log-density with respect to the parameters \\\mu\\ and \\\sigma\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma}
\tanh\left(\dfrac{y-\mu}{2\sigma}\right)\$\$ \$\$\dfrac{\partial
\ell}{\partial \sigma} = -\dfrac{1}{\sigma} \left\[ 1 -
\dfrac{y-\mu}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right)
\right\]\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of first derivatives.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
