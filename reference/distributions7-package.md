# distributions7: S7 Framework for Probability Distributions

Probability distributions as objects, built on the S7 object system.
Each distribution carries its density, distribution and quantile
functions, a random generator, its moments and the exact derivatives of
the log-likelihood up to the fourth order (and the fifth in closed form
for some families) with respect to the parameters, the response and the
parameters on the scale of their link functions. The package includes
univariate and multivariate families, wrappers for zero inflation, zero
adjustment, truncation, folding, transformation, fixed parameters and
reparametrization, maximum likelihood fitting with fit_distrib() and a
validation battery with check_distrib(). A distribution defined by its
density alone receives numerical methods for every other quantity. Part
of the 'statmodels7' toolkit.

## See also

Useful links:

- <https://statmodels7.github.io/distributions7/>

- <https://github.com/statmodels7/distributions7>

- Report bugs at <https://github.com/statmodels7/distributions7/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>
