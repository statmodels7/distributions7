# Package index

## Distributions

Each constructor returns a distribution object, taking the link function
used for each of its parameters.

- [`gaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
  : Gaussian Distribution Object (Standard Deviation Parameterization)
- [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
  : Cauchy Distribution Object
- [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
  : Logistic Distribution Object
- [`student_t_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
  : Student's t Distribution Object (Location-Scale Parameterization)
- [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
  : Laplace Distribution Object
- [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
  : Pseudo-Huber Distribution Object (Location-Scale Parameterization)
- [`skewnormal_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
  : Skew Normal Distribution Object
- [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
  : Skew t Distribution Object
- [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
  : Gumbel Distribution Object
- [`gamma_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
  : Gamma Distribution Object (Mean-Variance Parameterization)
- [`invgauss_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
  : Inverse-Gaussian Distribution Object (Mean-Dispersion
  Parameterization)
- [`lognormal_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
  : Lognormal Distribution Object (Log-Scale Parameterization)
- [`weibull_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
  : Weibull Distribution Object
- [`beta_distrib()`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
  : Beta Distribution Object (Mean-Precision Parameterization)
- [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
  : Bernoulli Distribution Object
- [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
  : Binomial Distribution Object
- [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
  : Poisson Distribution Object
- [`negbin_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
  : Negative Binomial Distribution Object (NB2)

## Probability functions

Density, distribution and quantile functions, and random generation.

- [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
  : Probability Density Function
- [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
  : Cumulative Distribution Function
- [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
  : Quantile Function
- [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  : Random Number Generator
- [`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
  : Atoms of a Distribution
- [`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
  : Generalized Ratio-of-Uniforms Sampling

## Derivatives

The score, the observed and expected information, and higher orders —
all available on the parameter scale or, through `scale = "link"`, with
respect to the unconstrained parameters.

- [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
  : Analytical Gradient
- [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
  : Analytical Hessian
- [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
  : Analytical Expected Hessian
- [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
  : Analytical Third-Order Derivatives
- [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
  : Analytical Fourth-Order Derivatives
- [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  : Gradient of the Log-Density with Respect to the Response
- [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
  : Second Derivative of the Log-Density with Respect to the Response
- [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  : Mixed Response-Parameter Derivatives of the Log-Density
- [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
  : Gradient of the Log Distribution Function
- [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
  : Second Derivatives of the Log Distribution Function
- [`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
  : Derivatives on the Link (Real) Scale
- [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  : Strategies for Expected Derivatives

## Estimation

- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  : Maximum-Likelihood Estimation
- [`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
  : Fisher Scoring, With Its Own Settings
- [`FisherScoring()`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md)
  : Fisher Scoring as an Object
- [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
  : A Starting Value Drawn from the Data
- [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  : S7 Class for Maximum-Likelihood Fits
- [`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
  : Numerically Validate a Distribution

## Moments

- [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  : Expected Value of a Function of a Random Variable
- [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  : Raw and Central Moments of a Distribution
- [`mean.distrib`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
  : Mean of a Distribution Object
- [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
  : Variance of a Distribution or Sample
- [`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
  : Standard Deviation of a Distribution or Sample
- [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
  : Skewness of a Distribution or Sample
- [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
  : Excess Kurtosis of a Distribution or Sample

## Building on a distribution

Wrappers that turn a distribution into another one: a mixture with a
point mass at zero, a restriction to an interval, the law of a
transformed variable, or the same law with some parameters held at known
values.

- [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  : Zero-Inflated Distribution Object (Discrete)
- [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
  : Zero-Adjusted Distribution Object
- [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
  : Truncated Distribution Object
- [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
  : Fix Parameters of a Distribution at Known Values
- [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
  : Apply a Variable Transformation to a Distribution Object
- [`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
  : S7 Class for Variable Transformers
- [`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md)
  : Logarithmic Transformation
- [`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md)
  : Affine (Location-Scale) Transformation
- [`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md)
  : Inverse Hyperbolic Sine Transformation
- [`bc_transform()`](https://statmodels7.github.io/distributions7/reference/bc_transform.md)
  : Box-Cox Transformation
- [`exp_transform()`](https://statmodels7.github.io/distributions7/reference/exp_transform.md)
  : Exponential Transformation
- [`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md)
  : Expit (Sigmoid) Transformation
- [`inverse_transform()`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md)
  : Reciprocal (Inverse) Transformation
- [`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md)
  : Logit Transformation
- [`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md)
  : Power Transformation
- [`softplus_transform()`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md)
  : Softplus Transformation
- [`sqrt_transform()`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md)
  : Square Root Transformation
- [`yj_transform()`](https://statmodels7.github.io/distributions7/reference/yj_transform.md)
  : Yeo-Johnson Transformation

## Several dimensions

Distributions whose observations are vectors. The response is a matrix
with one row per observation, and the parameters stay scalars: a mean
vector contributes one each, and a covariance contributes the free
values of the covstructs7 structure that parametrises it, so every
generic of the package indexes them as it always did.

- [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)
  : Construct a Multivariate Gaussian Distribution
- [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
  : Construct a Multivariate Student's t Distribution
- [`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
  [`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
  : The Mean Vector and Covariance a Parameter List Describes
- [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
  : A Marginal of a Multivariate Distribution
- [`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
  : Interpretable Quantities of a Multivariate Distribution
- [`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
  : Interpretable Estimates of a Multivariate Fit
- [`n_obs()`](https://statmodels7.github.io/distributions7/reference/n_obs.md)
  : How Many Observations a Response Holds

## Classes

The S7 classes; each page lists the methods that dispatch on it.

- [`distrib()`](https://statmodels7.github.io/distributions7/reference/distrib.md)
  : S7 Class for Probability Distributions
- [`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
  : S7 Class for Continuous Distributions
- [`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
  : S7 Class for Discrete Distributions
- [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  : S7 Class for Multivariate Distributions

## Numerical fallbacks

What a distribution gets for free when it implements only its density.
Rarely called directly, but useful as a reference for what is happening.

- [`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)
  : Numerical Gradient of the Log-Density
- [`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md)
  : Numerical Hessian of the Log-Density
- [`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md)
  : Numerical Third-Order Derivatives of the Log-Density
- [`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
  : Numerical Fourth-Order Derivatives of the Log-Density
- [`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)
  : Numerical Gradient of the Log-Density with Respect to the Response
- [`numerical_hess_y()`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md)
  : Numerical Second Derivative of the Log-Density with Respect to the
  Response
- [`numerical_cross_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)
  : Numerical Mixed Response-Parameter Derivatives
- [`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
  : Numerical Derivatives of the Distribution Function
- [`numerical_series()`](https://statmodels7.github.io/distributions7/reference/numerical_series.md)
  : Numerical Summation of Discrete Series

## Utilities

- [`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md)
  : Check Parameter Values Against Their Domains
- [`check_params_dim()`](https://statmodels7.github.io/distributions7/reference/check_params_dim.md)
  : Check Consistency of Parameter Dimensions
- [`expand_params()`](https://statmodels7.github.io/distributions7/reference/expand_params.md)
  : Expand Parameters to Common Length
- [`transpose_params()`](https://statmodels7.github.io/distributions7/reference/transpose_params.md)
  : Transpose and Simplify Parameter List Structure
- [`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
  : Generate Names for Higher-Order Derivative Components
- [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
  : Generate Names for Hessian Matrix Components
- [`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)
  : Per-Parameter Smoothness of the Log-Likelihood
- [`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
  : Generate Random Parameters
- [`distrib_generics`](https://statmodels7.github.io/distributions7/reference/distrib_generics.md)
  : Distribution Generics
- [`coef.distrib_fit`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md)
  : Extract Estimates from a Maximum-Likelihood Fit
- [`vcov.distrib_fit`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
  : Variance-Covariance Matrix of a Maximum-Likelihood Fit
- [`confint.distrib_fit`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
  : Confidence Intervals for a Maximum-Likelihood Fit
- [`logLik.distrib_fit`](https://statmodels7.github.io/distributions7/reference/logLik.distrib_fit.md)
  : Log-Likelihood of a Maximum-Likelihood Fit
- [`simulate.distrib_fit`](https://statmodels7.github.io/distributions7/reference/simulate.distrib_fit.md)
  : Simulate from a Fitted Distribution
- [`print.distrib_fit`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
  : Print Method for Maximum-Likelihood Fits
- [`plot.distrib_fit`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
  : Plot a Fitted Distribution against the Data

## Gaussian

- [`GaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/GaussianDistrib.md)
  : S7 Class for Gaussian Distribution
- [`MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  : Multivariate Gaussian Distribution
- [`distrib_cdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GaussianDistrib.md)
  : Gaussian Cumulative Distribution Function
- [`distrib_deriv3.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GaussianDistrib.md)
  : Gaussian Analytical Third-Order Derivatives
- [`distrib_deriv4.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GaussianDistrib.md)
  : Gaussian Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GaussianDistrib.md)
  : Gaussian Analytical Expected Hessian
- [`distrib_grad_cdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GaussianDistrib.md)
  : Gaussian Log-CDF Derivatives
- [`distrib_grad_y.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GaussianDistrib.md)
  : Gaussian Response Derivatives
- [`distrib_gradient.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GaussianDistrib.md)
  : Gaussian Analytical Gradient
- [`distrib_hess_cdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.GaussianDistrib.md)
  : Gaussian Log-CDF Second Derivatives
- [`distrib_hess_y.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.GaussianDistrib.md)
  : Gaussian Response Second Derivative
- [`distrib_hessian.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GaussianDistrib.md)
  : Gaussian Analytical Observed Hessian
- [`distrib_pdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GaussianDistrib.md)
  : Gaussian Probability Density Function
- [`distrib_quantile.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GaussianDistrib.md)
  : Gaussian Quantile Function
- [`distrib_rng.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GaussianDistrib.md)
  : Gaussian Random Number Generator

## Cauchy

- [`CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/CauchyDistrib.md)
  : S7 Class for Cauchy Distribution
- [`distrib_cdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.CauchyDistrib.md)
  : Cauchy Cumulative Distribution Function
- [`distrib_deriv3.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.CauchyDistrib.md)
  : Cauchy Analytical Third-Order Derivatives
- [`distrib_deriv4.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.CauchyDistrib.md)
  : Cauchy Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md)
  : Cauchy Analytical Expected Hessian
- [`distrib_grad_cdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.CauchyDistrib.md)
  : Cauchy Log-CDF Derivatives
- [`distrib_grad_y.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.CauchyDistrib.md)
  : Cauchy Response Derivatives
- [`distrib_gradient.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
  : Cauchy Analytical Gradient
- [`distrib_hess_cdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.CauchyDistrib.md)
  : Cauchy Log-CDF Second Derivatives
- [`distrib_hess_y.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.CauchyDistrib.md)
  : Cauchy Response Second Derivative
- [`distrib_hessian.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
  : Cauchy Analytical Observed Hessian
- [`distrib_pdf.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.CauchyDistrib.md)
  : Cauchy Probability Density Function
- [`distrib_quantile.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.CauchyDistrib.md)
  : Cauchy Quantile Function
- [`distrib_rng.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.CauchyDistrib.md)
  : Cauchy Random Number Generator

## Logistic

- [`LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/LogisticDistrib.md)
  : S7 Class for Logistic Distribution
- [`distrib_cdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LogisticDistrib.md)
  : Logistic Cumulative Distribution Function
- [`distrib_deriv3.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md)
  : Logistic Analytical Third-Order Derivatives
- [`distrib_deriv4.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LogisticDistrib.md)
  : Logistic Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md)
  : Logistic Analytical Expected Hessian
- [`distrib_grad_cdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LogisticDistrib.md)
  : Logistic Log-CDF Derivatives
- [`distrib_grad_y.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LogisticDistrib.md)
  : Logistic Response Derivatives
- [`distrib_gradient.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
  : Logistic Analytical Gradient
- [`distrib_hess_cdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LogisticDistrib.md)
  : Logistic Log-CDF Second Derivatives
- [`distrib_hess_y.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LogisticDistrib.md)
  : Logistic Response Second Derivative
- [`distrib_hessian.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
  : Logistic Analytical Observed Hessian
- [`distrib_pdf.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LogisticDistrib.md)
  : Logistic Probability Density Function
- [`distrib_quantile.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LogisticDistrib.md)
  : Logistic Quantile Function
- [`distrib_rng.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LogisticDistrib.md)
  : Logistic Random Number Generator

## Student’s t

- [`StudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/StudentTDistrib.md)
  : S7 Class for Student's t Distribution
- [`MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  : Multivariate Student's t Distribution
- [`distrib_cdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentTDistrib.md)
  : Student's t Cumulative Distribution Function
- [`distrib_deriv3.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentTDistrib.md)
  : Student's t Analytical Third-Order Derivatives
- [`distrib_deriv4.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentTDistrib.md)
  : Student's t Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentTDistrib.md)
  : Student's t Analytical Expected Hessian
- [`distrib_grad_cdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentTDistrib.md)
  : Student t Log-CDF Gradient
- [`distrib_grad_y.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.StudentTDistrib.md)
  : Student's t Response Derivatives
- [`distrib_gradient.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentTDistrib.md)
  : Student's t Analytical Gradient
- [`distrib_hess_y.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.StudentTDistrib.md)
  : Student's t Response Second Derivative
- [`distrib_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentTDistrib.md)
  : Student's t Analytical Observed Hessian
- [`distrib_pdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentTDistrib.md)
  : Student's t Probability Density Function
- [`distrib_quantile.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.StudentTDistrib.md)
  : Student's t Quantile Function
- [`distrib_rng.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.StudentTDistrib.md)
  : Student's t Random Number Generator

## Laplace

- [`LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/LaplaceDistrib.md)
  : S7 Class for Laplace Distribution
- [`distrib_cdf.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md)
  : Laplace Cumulative Distribution Function
- [`distrib_deriv3.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LaplaceDistrib.md)
  : Laplace Analytical Third-Order Derivatives
- [`distrib_deriv4.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LaplaceDistrib.md)
  : Laplace Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
  : Laplace Analytical Expected Hessian (Fisher Information)
- [`distrib_grad_cdf.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LaplaceDistrib.md)
  : Laplace Log-CDF Derivatives
- [`distrib_grad_y.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md)
  : Laplace Response Derivatives
- [`distrib_gradient.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md)
  : Laplace Analytical Gradient
- [`distrib_hess_cdf.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LaplaceDistrib.md)
  : Laplace Log-CDF Second Derivatives
- [`distrib_hess_y.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LaplaceDistrib.md)
  : Laplace Response Second Derivative
- [`distrib_hessian.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
  : Laplace Analytical Observed Hessian
- [`distrib_pdf.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LaplaceDistrib.md)
  : Laplace Probability Density Function
- [`distrib_quantile.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LaplaceDistrib.md)
  : Laplace Quantile Function
- [`distrib_rng.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LaplaceDistrib.md)
  : Laplace Random Number Generator

## Pseudo-Huber

- [`PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/PseudoHuberDistrib.md)
  : S7 Class for Pseudo-Huber Distribution
- [`distrib_cdf.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuberDistrib.md)
  : Pseudo-Huber Cumulative Distribution Function
- [`distrib_deriv3.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md)
  : Pseudo-Huber Analytical Third-Order Derivatives
- [`distrib_deriv4.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.PseudoHuberDistrib.md)
  : Pseudo-Huber Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
  : Pseudo-Huber Analytical Expected Hessian
- [`distrib_grad_cdf.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md)
  : Pseudo-Huber Log-CDF Gradient
- [`distrib_grad_y.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuberDistrib.md)
  : Pseudo-Huber Response Derivatives
- [`distrib_gradient.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)
  : Pseudo-Huber Analytical Gradient
- [`distrib_hess_y.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.PseudoHuberDistrib.md)
  : Pseudo-Huber Response Second Derivative
- [`distrib_hessian.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
  : Pseudo-Huber Analytical Observed Hessian
- [`distrib_pdf.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md)
  : Pseudo-Huber Probability Density Function
- [`distrib_quantile.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md)
  : Pseudo-Huber Quantile Function
- [`distrib_rng.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md)
  : Pseudo-Huber Random Number Generator

## Skew normal

- [`SkewNormalDistrib()`](https://statmodels7.github.io/distributions7/reference/SkewNormalDistrib.md)
  : S7 Class for the Skew Normal Distribution
- [`distrib_cdf.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.SkewNormalDistrib.md)
  : Skew Normal Cumulative Distribution Function
- [`distrib_deriv3.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormalDistrib.md)
  : Skew Normal Analytical Third-Order Derivatives
- [`distrib_deriv4.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewNormalDistrib.md)
  : Skew Normal Analytical Fourth-Order Derivatives
- [`distrib_grad_y.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewNormalDistrib.md)
  : Skew Normal Response Derivative
- [`distrib_gradient.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormalDistrib.md)
  : Skew Normal Analytical Gradient
- [`distrib_hess_y.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormalDistrib.md)
  : Skew Normal Response Second Derivative
- [`distrib_hessian.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormalDistrib.md)
  : Skew Normal Analytical Observed Hessian
- [`distrib_pdf.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormalDistrib.md)
  : Skew Normal Probability Density Function
- [`distrib_rng.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormalDistrib.md)
  : Skew Normal Random Number Generator

## Skew t

- [`SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/SkewTDistrib.md)
  : S7 Class for the Skew t Distribution
- [`distrib_deriv3.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md)
  : Skew t Third-Order Derivatives
- [`distrib_deriv4.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
  : Skew t Fourth-Order Derivatives
- [`distrib_grad_y.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewTDistrib.md)
  : Skew t Response Derivative
- [`distrib_gradient.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
  : Skew t Analytical Gradient
- [`distrib_hess_y.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewTDistrib.md)
  : Skew t Response Second Derivative
- [`distrib_hessian.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
  : Skew t Analytical Observed Hessian
- [`distrib_pdf.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewTDistrib.md)
  : Skew t Probability Density Function
- [`distrib_rng.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewTDistrib.md)
  : Skew t Random Number Generator

## Gumbel

- [`GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/GumbelDistrib.md)
  : S7 Class for the Gumbel Distribution
- [`distrib_cdf.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GumbelDistrib.md)
  : Gumbel Cumulative Distribution Function
- [`distrib_deriv3.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md)
  : Gumbel Analytical Third-Order Derivatives
- [`distrib_deriv4.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GumbelDistrib.md)
  : Gumbel Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md)
  : Gumbel Analytical Expected Hessian
- [`distrib_grad_y.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GumbelDistrib.md)
  : Gumbel Response Derivative
- [`distrib_gradient.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md)
  : Gumbel Analytical Gradient
- [`distrib_hess_y.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.GumbelDistrib.md)
  : Gumbel Response Second Derivative
- [`distrib_hessian.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md)
  : Gumbel Analytical Observed Hessian
- [`distrib_pdf.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GumbelDistrib.md)
  : Gumbel Probability Density Function
- [`distrib_quantile.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GumbelDistrib.md)
  : Gumbel Quantile Function
- [`distrib_rng.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GumbelDistrib.md)
  : Gumbel Random Number Generator

## Weibull

- [`WeibullDistrib()`](https://statmodels7.github.io/distributions7/reference/WeibullDistrib.md)
  : S7 Class for the Weibull Distribution
- [`distrib_cdf.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.WeibullDistrib.md)
  : Weibull Cumulative Distribution Function
- [`distrib_deriv3.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.WeibullDistrib.md)
  : Weibull Analytical Third-Order Derivatives
- [`distrib_deriv4.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.WeibullDistrib.md)
  : Weibull Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.WeibullDistrib.md)
  : Weibull Analytical Expected Hessian
- [`distrib_grad_y.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.WeibullDistrib.md)
  : Weibull Response Derivative
- [`distrib_gradient.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.WeibullDistrib.md)
  : Weibull Analytical Gradient
- [`distrib_hess_y.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.WeibullDistrib.md)
  : Weibull Response Second Derivative
- [`distrib_hessian.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.WeibullDistrib.md)
  : Weibull Analytical Observed Hessian
- [`distrib_pdf.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.WeibullDistrib.md)
  : Weibull Probability Density Function
- [`distrib_quantile.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.WeibullDistrib.md)
  : Weibull Quantile Function
- [`distrib_rng.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.WeibullDistrib.md)
  : Weibull Random Number Generator

## Gamma

- [`GammaDistrib()`](https://statmodels7.github.io/distributions7/reference/GammaDistrib.md)
  : S7 Class for Gamma Distribution
- [`distrib_cdf.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GammaDistrib.md)
  : Gamma Cumulative Distribution Function
- [`distrib_deriv3.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GammaDistrib.md)
  : Gamma Analytical Third-Order Derivatives
- [`distrib_deriv4.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GammaDistrib.md)
  : Gamma Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GammaDistrib.md)
  : Gamma Analytical Expected Hessian
- [`distrib_grad_y.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GammaDistrib.md)
  : Gamma Response Derivatives
- [`distrib_gradient.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GammaDistrib.md)
  : Gamma Analytical Gradient
- [`distrib_hess_y.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.GammaDistrib.md)
  : Gamma Response Second Derivative
- [`distrib_hessian.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GammaDistrib.md)
  : Gamma Analytical Observed Hessian
- [`distrib_pdf.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GammaDistrib.md)
  : Gamma Probability Density Function
- [`distrib_quantile.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GammaDistrib.md)
  : Gamma Quantile Function
- [`distrib_rng.GammaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GammaDistrib.md)
  : Gamma Random Number Generator

## Inverse Gaussian

- [`InvGaussDistrib()`](https://statmodels7.github.io/distributions7/reference/InvGaussDistrib.md)
  : S7 Class for Inverse-Gaussian Distribution
- [`distrib_cdf.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGaussDistrib.md)
  : Inverse-Gaussian Cumulative Distribution Function
- [`distrib_deriv3.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGaussDistrib.md)
  : Inverse-Gaussian Analytical Third-Order Derivatives
- [`distrib_deriv4.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGaussDistrib.md)
  : Inverse-Gaussian Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGaussDistrib.md)
  : Inverse-Gaussian Analytical Expected Hessian
- [`distrib_grad_cdf.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGaussDistrib.md)
  : Inverse Gaussian Log-CDF Gradient
- [`distrib_grad_y.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGaussDistrib.md)
  : Inverse-Gaussian Response Derivatives
- [`distrib_gradient.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGaussDistrib.md)
  : Inverse-Gaussian Analytical Gradient
- [`distrib_hess_y.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.InvGaussDistrib.md)
  : Inverse-Gaussian Response Second Derivative
- [`distrib_hessian.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGaussDistrib.md)
  : Inverse-Gaussian Analytical Observed Hessian
- [`distrib_pdf.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGaussDistrib.md)
  : Inverse-Gaussian Probability Density Function
- [`distrib_quantile.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGaussDistrib.md)
  : Inverse-Gaussian Quantile Function
- [`distrib_rng.InvGaussDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.InvGaussDistrib.md)
  : Inverse-Gaussian Random Number Generator

## Lognormal

- [`LognormalDistrib()`](https://statmodels7.github.io/distributions7/reference/LognormalDistrib.md)
  : S7 Class for Lognormal Distribution
- [`distrib_cdf.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LognormalDistrib.md)
  : Lognormal Cumulative Distribution Function
- [`distrib_deriv3.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LognormalDistrib.md)
  : Lognormal Analytical Third-Order Derivatives
- [`distrib_deriv4.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LognormalDistrib.md)
  : Lognormal Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LognormalDistrib.md)
  : Lognormal Analytical Expected Hessian
- [`distrib_grad_cdf.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LognormalDistrib.md)
  : Lognormal Log-CDF Gradient
- [`distrib_grad_y.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LognormalDistrib.md)
  : Lognormal Response Derivatives
- [`distrib_gradient.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LognormalDistrib.md)
  : Lognormal Analytical Gradient
- [`distrib_hess_y.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LognormalDistrib.md)
  : Lognormal Response Second Derivative
- [`distrib_hessian.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LognormalDistrib.md)
  : Lognormal Analytical Observed Hessian
- [`distrib_pdf.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LognormalDistrib.md)
  : Lognormal Probability Density Function
- [`distrib_quantile.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LognormalDistrib.md)
  : Lognormal Quantile Function
- [`distrib_rng.LognormalDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LognormalDistrib.md)
  : Lognormal Random Number Generator

## Beta

- [`BetaDistrib()`](https://statmodels7.github.io/distributions7/reference/BetaDistrib.md)
  : S7 Class for Beta Distribution
- [`distrib_cdf.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaDistrib.md)
  : Beta Cumulative Distribution Function
- [`distrib_deriv3.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaDistrib.md)
  : Beta Analytical Third-Order Derivatives
- [`distrib_deriv4.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaDistrib.md)
  : Beta Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaDistrib.md)
  : Beta Analytical Expected Hessian
- [`distrib_grad_y.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.BetaDistrib.md)
  : Beta Response Derivatives
- [`distrib_gradient.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaDistrib.md)
  : Beta Analytical Gradient
- [`distrib_hess_y.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.BetaDistrib.md)
  : Beta Response Second Derivative
- [`distrib_hessian.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaDistrib.md)
  : Beta Analytical Observed Hessian
- [`distrib_pdf.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaDistrib.md)
  : Beta Probability Density Function
- [`distrib_quantile.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BetaDistrib.md)
  : Beta Quantile Function
- [`distrib_rng.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaDistrib.md)
  : Beta Random Number Generator

## Bernoulli

- [`BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/BernoulliDistrib.md)
  : S7 Class for Bernoulli Distribution
- [`distrib_cdf.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BernoulliDistrib.md)
  : Bernoulli Cumulative Distribution Function
- [`distrib_deriv3.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BernoulliDistrib.md)
  : Bernoulli Analytical Third-Order Derivatives
- [`distrib_deriv4.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BernoulliDistrib.md)
  : Bernoulli Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BernoulliDistrib.md)
  : Bernoulli Analytical Expected Hessian
- [`distrib_grad_cdf.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.BernoulliDistrib.md)
  : Bernoulli Log-CDF Gradient
- [`distrib_gradient.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BernoulliDistrib.md)
  : Bernoulli Analytical Gradient
- [`distrib_hessian.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BernoulliDistrib.md)
  : Bernoulli Analytical Observed Hessian
- [`distrib_pdf.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BernoulliDistrib.md)
  : Bernoulli Probability Mass Function
- [`distrib_quantile.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BernoulliDistrib.md)
  : Bernoulli Quantile Function
- [`distrib_rng.BernoulliDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BernoulliDistrib.md)
  : Bernoulli Random Number Generator

## Binomial

- [`BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/BinomialDistrib.md)
  : S7 Class for Binomial Distribution
- [`distrib_cdf.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md)
  : Binomial Cumulative Distribution Function
- [`distrib_deriv3.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BinomialDistrib.md)
  : Binomial Analytical Third-Order Derivatives
- [`distrib_deriv4.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BinomialDistrib.md)
  : Binomial Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md)
  : Binomial Analytical Expected Hessian
- [`distrib_grad_cdf.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.BinomialDistrib.md)
  : Binomial Log-CDF Gradient
- [`distrib_gradient.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BinomialDistrib.md)
  : Binomial Analytical Gradient
- [`distrib_hessian.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BinomialDistrib.md)
  : Binomial Analytical Observed Hessian
- [`distrib_pdf.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BinomialDistrib.md)
  : Binomial Probability Mass Function
- [`distrib_quantile.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BinomialDistrib.md)
  : Binomial Quantile Function
- [`distrib_rng.BinomialDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BinomialDistrib.md)
  : Binomial Random Number Generator

## Poisson

- [`PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/PoissonDistrib.md)
  : S7 Class for Poisson Distribution
- [`distrib_cdf.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md)
  : Poisson Cumulative Distribution Function
- [`distrib_deriv3.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PoissonDistrib.md)
  : Poisson Analytical Third-Order Derivatives
- [`distrib_deriv4.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.PoissonDistrib.md)
  : Poisson Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
  : Poisson Analytical Expected Hessian
- [`distrib_grad_cdf.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PoissonDistrib.md)
  : Poisson Log-CDF Gradient
- [`distrib_gradient.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PoissonDistrib.md)
  : Poisson Analytical Gradient
- [`distrib_hessian.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PoissonDistrib.md)
  : Poisson Analytical Observed Hessian
- [`distrib_pdf.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PoissonDistrib.md)
  : Poisson Probability Mass Function
- [`distrib_quantile.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PoissonDistrib.md)
  : Poisson Quantile Function
- [`distrib_rng.PoissonDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PoissonDistrib.md)
  : Poisson Random Number Generator

## Negative binomial

- [`NegBinDistrib()`](https://statmodels7.github.io/distributions7/reference/NegBinDistrib.md)
  : S7 Class for Negative Binomial Distribution (NB2)
- [`distrib_cdf.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBinDistrib.md)
  : Negative Binomial Cumulative Distribution Function
- [`distrib_deriv3.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBinDistrib.md)
  : Negative Binomial Analytical Third-Order Derivatives
- [`distrib_deriv4.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBinDistrib.md)
  : Negative Binomial Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBinDistrib.md)
  : Negative Binomial Analytical Expected Hessian
- [`distrib_grad_cdf.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.NegBinDistrib.md)
  : Negative Binomial Log-CDF Gradient
- [`distrib_gradient.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBinDistrib.md)
  : Negative Binomial Analytical Gradient
- [`distrib_hessian.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBinDistrib.md)
  : Negative Binomial Analytical Observed Hessian
- [`distrib_pdf.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBinDistrib.md)
  : Negative Binomial Probability Mass Function
- [`distrib_quantile.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBinDistrib.md)
  : Negative Binomial Quantile Function
- [`distrib_rng.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBinDistrib.md)
  : Negative Binomial Random Number Generator

## Multivariate gaussian

- [`MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  : Multivariate Gaussian Distribution

## Multivariate Student t

- [`MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  : Multivariate Student's t Distribution

## Wrapper classes

The methods of the zero-inflated, zero-adjusted, truncated, transformed
and fixed-parameter wrappers.

- [`ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md)
  : S7 Class for Zero-Inflated Distributions
- [`distrib_cdf.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroInflatedDistrib.md)
  : Zero-Inflated Cumulative Distribution Function
- [`distrib_deriv3.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroInflatedDistrib.md)
  : Zero-Inflated Third Derivatives
- [`distrib_deriv4.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroInflatedDistrib.md)
  : Zero-Inflated Fourth Derivatives
- [`distrib_expected_hessian.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroInflatedDistrib.md)
  : Zero-Inflated Analytical Expected Hessian
- [`distrib_gradient.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroInflatedDistrib.md)
  : Zero-Inflated Analytical Gradient
- [`distrib_hessian.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroInflatedDistrib.md)
  : Zero-Inflated Analytical Observed Hessian
- [`distrib_pdf.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroInflatedDistrib.md)
  : Zero-Inflated Probability Mass Function
- [`distrib_quantile.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroInflatedDistrib.md)
  : Zero-Inflated Quantile Function
- [`distrib_rng.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroInflatedDistrib.md)
  : Zero-Inflated Random Number Generator
- [`ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedDiscreteDistrib.md)
  : S7 Class for Zero-Adjusted Discrete (Hurdle) Distributions
- [`distrib_cdf.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Cumulative Distribution Function
- [`distrib_deriv3.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroAdjustedDiscreteDistrib.md)
  : Hurdle Third Derivatives
- [`distrib_deriv4.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroAdjustedDiscreteDistrib.md)
  : Hurdle Fourth Derivatives
- [`distrib_expected_hessian.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Analytical Expected Hessian
- [`distrib_gradient.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Analytical Gradient
- [`distrib_hessian.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Analytical Observed Hessian
- [`distrib_pdf.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Probability Mass Function
- [`distrib_quantile.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Quantile Function
- [`distrib_rng.ZeroAdjustedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroAdjustedDiscreteDistrib.md)
  : Zero-Adjusted Discrete Random Number Generator
- [`ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/ZeroAdjustedContinuousDistrib.md)
  : S7 Class for Zero-Adjusted Continuous Distributions
- [`distrib_atoms.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
  : Atoms of a Zero-Adjusted Continuous Distribution
- [`distrib_cdf.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Cumulative Distribution Function
- [`distrib_deriv3.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Third Derivatives
- [`distrib_deriv4.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Fourth Derivatives
- [`distrib_expected_hessian.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Analytical Expected Hessian
- [`distrib_grad_y.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Response Gradient
- [`distrib_gradient.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Analytical Gradient
- [`distrib_hess_y.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Response Hessian
- [`distrib_hessian.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Analytical Observed Hessian
- [`distrib_pdf.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Probability Density Function
- [`distrib_quantile.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Quantile Function
- [`distrib_rng.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroAdjustedContinuousDistrib.md)
  : Zero-Adjusted Continuous Random Number Generator
- [`TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/TruncatedDiscreteDistrib.md)
  : S7 Class for Truncated Discrete Distributions
- [`distrib_cdf.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md)
  : Truncated Cumulative Distribution Function (Discrete)
- [`distrib_deriv3.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TruncatedDiscreteDistrib.md)
  : Truncated Third Derivatives (Discrete)
- [`distrib_deriv4.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TruncatedDiscreteDistrib.md)
  : Truncated Fourth Derivatives (Discrete)
- [`distrib_expected_hessian.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedDiscreteDistrib.md)
  : Truncated Analytical Expected Hessian (Discrete)
- [`distrib_gradient.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedDiscreteDistrib.md)
  : Truncated Analytical Gradient (Discrete)
- [`distrib_hessian.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedDiscreteDistrib.md)
  : Truncated Analytical Observed Hessian (Discrete)
- [`distrib_pdf.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedDiscreteDistrib.md)
  : Truncated Probability Mass Function
- [`distrib_quantile.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedDiscreteDistrib.md)
  : Truncated Quantile Function (Discrete)
- [`distrib_rng.TruncatedDiscreteDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedDiscreteDistrib.md)
  : Truncated Random Number Generator (Discrete)
- [`TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/TruncatedContinuousDistrib.md)
  : S7 Class for Truncated Continuous Distributions
- [`distrib_atoms.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md)
  : Atoms of a Truncated Continuous Distribution
- [`distrib_cdf.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md)
  : Truncated Cumulative Distribution Function (Continuous)
- [`distrib_deriv3.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TruncatedContinuousDistrib.md)
  : Truncated Third Derivatives (Continuous)
- [`distrib_deriv4.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TruncatedContinuousDistrib.md)
  : Truncated Fourth Derivatives (Continuous)
- [`distrib_expected_hessian.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedContinuousDistrib.md)
  : Truncated Analytical Expected Hessian (Continuous)
- [`distrib_grad_y.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.TruncatedContinuousDistrib.md)
  : Truncated Continuous Response Gradient
- [`distrib_gradient.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedContinuousDistrib.md)
  : Truncated Analytical Gradient (Continuous)
- [`distrib_hess_y.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.TruncatedContinuousDistrib.md)
  : Truncated Continuous Response Hessian
- [`distrib_hessian.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedContinuousDistrib.md)
  : Truncated Analytical Observed Hessian (Continuous)
- [`distrib_pdf.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedContinuousDistrib.md)
  : Truncated Probability Density Function
- [`distrib_quantile.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedContinuousDistrib.md)
  : Truncated Quantile Function (Continuous)
- [`distrib_rng.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedContinuousDistrib.md)
  : Truncated Random Number Generator (Continuous)
- [`TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md)
  : S7 Class for Transformed Distributions
- [`distrib_cdf.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TransformedDistrib.md)
  : Transformed Cumulative Distribution Function
- [`distrib_deriv3.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TransformedDistrib.md)
  : Transformed Third Derivatives
- [`distrib_deriv4.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TransformedDistrib.md)
  : Transformed Fourth Derivatives
- [`distrib_expected_hessian.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TransformedDistrib.md)
  : Transformed Analytical Expected Hessian
- [`distrib_gradient.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TransformedDistrib.md)
  : Transformed Analytical Gradient
- [`distrib_hessian.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TransformedDistrib.md)
  : Transformed Analytical Observed Hessian
- [`distrib_pdf.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TransformedDistrib.md)
  : Transformed Probability Density Function
- [`distrib_quantile.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TransformedDistrib.md)
  : Transformed Quantile Function
- [`distrib_rng.TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TransformedDistrib.md)
  : Transformed Random Number Generator
- [`FixedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
  : S7 Class for Distributions With Fixed Parameters (Continuous)
- [`FixedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/FixedDiscreteDistrib.md)
  : S7 Class for Distributions With Fixed Parameters (Discrete)

## Base-class defaults

The methods registered on the base classes, which every distribution
inherits unless it registers something more specific.

- [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
  : Bernoulli Distribution Object
- [`beta_distrib()`](https://statmodels7.github.io/distributions7/reference/beta_distrib.md)
  : Beta Distribution Object (Mean-Precision Parameterization)
- [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
  : Binomial Distribution Object
- [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
  : Cauchy Distribution Object
- [`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
  : Numerically Validate a Distribution
- [`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
  : S7 Class for Continuous Distributions
- [`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
  : S7 Class for Discrete Distributions
- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  : Maximum-Likelihood Estimation
- [`gamma_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
  : Gamma Distribution Object (Mean-Variance Parameterization)
- [`gaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
  : Gaussian Distribution Object (Standard Deviation Parameterization)
- [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
  : Gumbel Distribution Object
- [`invgauss_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
  : Inverse-Gaussian Distribution Object (Mean-Dispersion
  Parameterization)
- [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
  : Laplace Distribution Object
- [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
  : Logistic Distribution Object
- [`lognormal_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
  : Lognormal Distribution Object (Log-Scale Parameterization)
- [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  : S7 Class for Multivariate Distributions
- [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)
  : Construct a Multivariate Gaussian Distribution
- [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
  : Construct a Multivariate Student's t Distribution
- [`negbin_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
  : Negative Binomial Distribution Object (NB2)
- [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
  : Poisson Distribution Object
- [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
  : Pseudo-Huber Distribution Object (Location-Scale Parameterization)
- [`skewnormal_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
  : Skew Normal Distribution Object
- [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
  : Skew t Distribution Object
- [`student_t_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
  : Student's t Distribution Object (Location-Scale Parameterization)
- [`weibull_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
  : Weibull Distribution Object

## Internals

The machinery the exported functions are built from. None of it is
exported, and none of it is needed to use the package — it is documented
so that the derivations can be followed from the code that implements
them: the set partitions behind the Bartlett identities and the wrapper
derivatives, the Bell polynomials behind the link scale, the
ratio-of-uniforms sampler and the transforms that make it work on a
divergent density.

- [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md)
  : Align Parameters to the Distribution's Parameter Order

- [`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md)
  : Coerce a Multivariate Response to a Matrix

- [`assemble_deriv()`](https://statmodels7.github.io/distributions7/reference/assemble_deriv.md)
  : Assemble One Order of a Wrapper's Derivatives

- [`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
  : Complete Bell Polynomial in the Parent's Log-Derivatives

- [`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md)
  : Partial Bell Polynomials for Orders up to Four

- [`canon_key()`](https://statmodels7.github.io/distributions7/reference/canon_key.md)
  : Canonical Component Name of a Block

- [`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
  : Put CDF Derivatives on the Requested Tail and Scale

- [`check_bounds_fast()`](https://statmodels7.github.io/distributions7/reference/check_bounds_fast.md)
  : Check Parameter Domains, Taking the Properties as Arguments

- [`check_derivative_args()`](https://statmodels7.github.io/distributions7/reference/check_derivative_args.md)
  : Shared Argument Handling for the Derivative Generics

- [`check_distrib_mv()`](https://statmodels7.github.io/distributions7/reference/check_distrib_mv.md)
  : Validate a Multivariate Distribution

- [`check_not_stacked()`](https://statmodels7.github.io/distributions7/reference/check_not_stacked.md)
  : Refuse to Stack Two Zero Parameters

- [`check_support_is_rich_enough()`](https://statmodels7.github.io/distributions7/reference/check_support_is_rich_enough.md)
  : Refuse a Model With More Parameters Than the Support Can Distinguish

- [`check_truncation_points()`](https://statmodels7.github.io/distributions7/reference/check_truncation_points.md)
  : Validate the Truncation Endpoints

- [`deriv_index_list()`](https://statmodels7.github.io/distributions7/reference/deriv_index_list.md)
  : Index Tuples Matching the Package's Component Naming

- [`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
  : Index Tuples Behind the Higher-Order Derivative Names

- [`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
  : Cumulative Probability Table for a Discrete Distribution

- [`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
  : CDF Derivatives of a Discrete Distribution

- [`distrib_atoms.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.distrib.md)
  : Default Atoms: None

- [`distrib_cdf.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md)
  : Default Numerical CDF for Continuous Distributions

- [`distrib_cdf.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md)
  : Default Numerical CDF for Discrete Distributions

- [`distrib_cdf.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.multivariate_distrib.md)
  : No Distribution Function in Several Dimensions

- [`distrib_cross_y.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.GaussianDistrib.md)
  : Gaussian Mixed Derivatives

- [`distrib_cross_y.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.StudentTDistrib.md)
  : Student's t Mixed Derivatives

- [`distrib_cross_y.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.TruncatedContinuousDistrib.md)
  : Mixed Derivatives of a Truncated Distribution

- [`distrib_cross_y.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.continuous_distrib.md)
  : Default Mixed Derivatives for Continuous Distributions

- [`distrib_deriv3.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.distrib.md)
  :

  Default Third-Order Derivatives for `distrib` Objects

- [`distrib_deriv4.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.distrib.md)
  :

  Default Fourth-Order Derivatives for `distrib` Objects

- [`distrib_deriv_component()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv_component.md)
  : One Component of the Parent's Derivative

- [`distrib_expected_hessian.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvGaussianDistrib.md)
  : Multivariate Gaussian Expected Hessian

- [`distrib_expected_hessian.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.distrib.md)
  :

  Default Expected Hessian for `distrib` Objects

- [`distrib_expected_hessian.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.multivariate_distrib.md)
  : Expected Information of a Multivariate Distribution

- [`distrib_grad_cdf.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.continuous_distrib.md)
  : Default Log-CDF Gradient for Continuous Distributions

- [`distrib_grad_cdf.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md)
  : Log-CDF Gradient for Discrete Distributions

- [`distrib_grad_y.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvGaussianDistrib.md)
  : Multivariate Gaussian Response Gradient

- [`distrib_grad_y.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.MvStudentTDistrib.md)
  : Multivariate Student t Response Gradient

- [`distrib_grad_y.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.continuous_distrib.md)
  : Default Response Gradient for Continuous Distributions

- [`distrib_grad_y.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.multivariate_distrib.md)
  [`distrib_hess_y.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.multivariate_distrib.md)
  [`distrib_cross_y.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.multivariate_distrib.md)
  : Response Derivatives of a Multivariate Distribution

- [`distrib_gradient.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
  : Multivariate Gaussian Score

- [`distrib_gradient.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
  : Multivariate Student t Score

- [`distrib_gradient.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.distrib.md)
  :

  Default Numerical Gradient for `distrib` Objects

- [`distrib_hess_cdf.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.continuous_distrib.md)
  : Default Log-CDF Hessian for Continuous Distributions

- [`distrib_hess_cdf.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.discrete_distrib.md)
  : Log-CDF Hessian for Discrete Distributions

- [`distrib_hess_y.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvGaussianDistrib.md)
  : Multivariate Gaussian Response Hessian

- [`distrib_hess_y.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md)
  : Multivariate Student t Response Hessian

- [`distrib_hess_y.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.continuous_distrib.md)
  : Default Response Hessian for Continuous Distributions

- [`distrib_hessian.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
  : Multivariate Gaussian Observed Hessian

- [`distrib_hessian.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md)
  : Multivariate Student t Observed Hessian

- [`distrib_hessian.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.distrib.md)
  :

  Default Numerical Hessian for `distrib` Objects

- [`distrib_pdf.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvGaussianDistrib.md)
  : Multivariate Gaussian Density

- [`distrib_pdf.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvStudentTDistrib.md)
  : Multivariate Student t Density

- [`distrib_quantile.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md)
  : Default Numerical Quantile Function for Continuous Distributions

- [`distrib_quantile.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md)
  : Default Numerical Quantile Function for Discrete Distributions

- [`distrib_quantile.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.multivariate_distrib.md)
  : No Quantile Function in Several Dimensions

- [`distrib_rng.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvGaussianDistrib.md)
  : Multivariate Gaussian Generator

- [`distrib_rng.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MvStudentTDistrib.md)
  : Multivariate Student t Generator

- [`distrib_rng.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.continuous_distrib.md)
  : Default Numerical RNG for Continuous Distributions

- [`distrib_rng.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.discrete_distrib.md)
  : Default Numerical RNG for Discrete Distributions

- [`distrib_start.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvGaussianDistrib.md)
  : The Maximum Likelihood Estimate as a Starting Value

- [`distrib_start.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_start.MvStudentTDistrib.md)
  : The Gaussian Estimate as a Starting Value for a t

- [`distrib_start.distrib`](https://statmodels7.github.io/distributions7/reference/distrib_start.distrib.md)
  : Random Starting Values

- [`expectation.TruncatedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/expectation.TruncatedContinuousDistrib.md)
  : Expectation for Truncated Continuous Distributions

- [`expectation.ZeroAdjustedContinuousDistrib`](https://statmodels7.github.io/distributions7/reference/expectation.ZeroAdjustedContinuousDistrib.md)
  : Expectation for Zero-Adjusted Continuous Distributions

- [`expectation.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/expectation.continuous_distrib.md)
  : Expectation of a Continuous Distribution

- [`expectation.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/expectation.discrete_distrib.md)
  : Expectation of a Discrete Distribution

- [`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
  : Expected Derivatives by the Bartlett Identity

- [`expected_by_integrate()`](https://statmodels7.github.io/distributions7/reference/expected_by_integrate.md)
  : Expected Derivatives by Numerical Integration

- [`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
  : Expected Derivatives by Monte Carlo

- [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  : Dispatch an Expected-Derivative Strategy

- [`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
  : A Five-Point First Derivative

- [`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md)
  : A Five-Point Fourth Derivative

- [`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md)
  : A Five-Point Second Derivative

- [`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md)
  : A Five-Point Third Derivative

- [`fd_is_reliable()`](https://statmodels7.github.io/distributions7/reference/fd_is_reliable.md)
  : Which Observations the Finite-Difference Reference Can Be Trusted At

- [`fd_second()`](https://statmodels7.github.io/distributions7/reference/fd_second.md)
  : A Second Derivative From One Stencil

- [`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
  : Finite-Difference Steps That Respect a Parameter's Domain

- [`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
  : Finite-Difference Steps That Respect the Support

- [`find_lp_anchor()`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md)
  : Locate a High-Density Point of a Bare Log-Density

- [`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md)
  : Locate an Interior High-Density Point

- [`fit_dtheta_deta()`](https://statmodels7.github.io/distributions7/reference/fit_dtheta_deta.md)
  : Jacobian of the Inverse Link at the Estimate

- [`fit_eta_from_theta()`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md)
  : Map Parameters to the Link Scale

- [`fit_format_elapsed()`](https://statmodels7.github.io/distributions7/reference/fit_format_elapsed.md)
  : Render a Duration With a Unit Matched to Its Size

- [`fit_hess_matrix()`](https://statmodels7.github.io/distributions7/reference/fit_hess_matrix.md)
  : Summed Hessian on the Link Scale, as a Matrix

- [`fit_loglik()`](https://statmodels7.github.io/distributions7/reference/fit_loglik.md)
  : Total Log-Likelihood

- [`fit_score()`](https://statmodels7.github.io/distributions7/reference/fit_score.md)
  : Summed Score on the Link Scale

- [`fit_theta_from_eta()`](https://statmodels7.github.io/distributions7/reference/fit_theta_from_eta.md)
  : Map the Link Scale Back to Parameters

- [`fixed_full_theta()`](https://statmodels7.github.io/distributions7/reference/fixed_full_theta.md)
  : Splice the Fixed Values Back Into a Full Parameter List

- [`generate_random_theta.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.MvGaussianDistrib.md)
  : Random Parameters for a Multivariate Gaussian

- [`generate_random_theta.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.MvStudentTDistrib.md)
  : Random Parameters for a Multivariate Student t

- [`generate_random_theta.distrib`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.distrib.md)
  :

  Generate Random Parameters for `distrib` Objects

- [`grou_core()`](https://statmodels7.github.io/distributions7/reference/grou_core.md)
  : The Generalized Ratio-of-Uniforms Sampler

- [`grou_two_sided()`](https://statmodels7.github.io/distributions7/reference/grou_two_sided.md)
  : Ratio-of-Uniforms for a Density Diverging at Both Edges

- [`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md)
  : Does This Distribution Have a Real Quantile Method?

- [`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
  : Can the Parent Supply Exact CDF Derivatives?

- [`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
  : Does This Distribution Compute Its Expected Information Exactly?

- [`hess_pairs()`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
  : Invert the Hessian Component Names

- [`index_partitions()`](https://statmodels7.github.io/distributions7/reference/index_partitions.md)
  : Set Partitions of a Multi-Index

- [`inverse_link_derivs()`](https://statmodels7.github.io/distributions7/reference/inverse_link_derivs.md)
  : Inverse-Link Derivatives for Every Parameter

- [`is_class()`](https://statmodels7.github.io/distributions7/reference/is_class.md)
  : Is an S7 Class the Given Base Class?

- [`is_fixed()`](https://statmodels7.github.io/distributions7/reference/is_fixed.md)
  : Is This a Fixed-Parameter Wrapper?

- [`is_truncated()`](https://statmodels7.github.io/distributions7/reference/is_truncated.md)
  : Is This Distribution Already Truncated?

- [`is_zero_wrapper()`](https://statmodels7.github.io/distributions7/reference/is_zero_wrapper.md)
  : Does This Distribution Already Model a Probability of Zero?

- [`kurtosis.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.GumbelDistrib.md)
  : Kurtosis of the Gumbel Distribution

- [`kurtosis.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md)
  : Kurtosis of the Laplace Distribution

- [`kurtosis.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBinDistrib.md)
  : Kurtosis of the Negative Binomial Distribution

- [`kurtosis.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.PseudoHuberDistrib.md)
  : Kurtosis of the Pseudo-Huber Distribution

- [`kurtosis.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormalDistrib.md)
  : Kurtosis of the Skew Normal Distribution

- [`kurtosis.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewTDistrib.md)
  : Kurtosis of the Skew t Distribution

- [`kurtosis.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.WeibullDistrib.md)
  : Kurtosis of the Weibull Distribution

- [`kurtosis.distrib`](https://statmodels7.github.io/distributions7/reference/kurtosis.distrib.md)
  : Kurtosis of a Distribution

- [`kurtosis.numeric`](https://statmodels7.github.io/distributions7/reference/kurtosis.numeric.md)
  : Sample Kurtosis

- [`link_scale_lower_orders()`](https://statmodels7.github.io/distributions7/reference/link_scale_lower_orders.md)
  : Lower-Order Parameter-Scale Derivatives for the Chain Rule

- [`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
  : CDF Derivatives of a Location-Scale Family

- [`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
  : Location-Scale CDF Gradient

- [`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
  : Location-Scale CDF Hessian

- [`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
  : Derivatives of a Logarithm From the Ratios Alone

- [`log_pow_deriv()`](https://statmodels7.github.io/distributions7/reference/log_pow_deriv.md)
  : Derivatives of log(p) and log(1 - p)

- [`lp_edge_divergence()`](https://statmodels7.github.io/distributions7/reference/lp_edge_divergence.md)
  : Detect and Measure a Divergence at the Edges of the Support

- [`mean.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/mean.GumbelDistrib.md)
  : Mean of the Gumbel Distribution

- [`mean.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md)
  : Mean of the Laplace Distribution

- [`mean.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/mean.MvGaussianDistrib.md)
  : Mean of a Multivariate Gaussian

- [`mean.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/mean.MvStudentTDistrib.md)
  : Mean of a Multivariate Student t

- [`mean.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/mean.NegBinDistrib.md)
  : Mean of the Negative Binomial Distribution

- [`mean.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/mean.PseudoHuberDistrib.md)
  : Mean of the Pseudo-Huber Distribution

- [`mean.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormalDistrib.md)
  : Mean of the Skew Normal Distribution

- [`mean.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/mean.SkewTDistrib.md)
  : Mean of the Skew t Distribution

- [`mean.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/mean.WeibullDistrib.md)
  : Mean of the Weibull Distribution

- [`mean.distrib`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md)
  : Mean of a Distribution Object

- [`memo_ratio()`](https://statmodels7.github.io/distributions7/reference/memo_ratio.md)
  : Memoise a Ratio Function on Its Block

- [`mills_ratio()`](https://statmodels7.github.io/distributions7/reference/mills_ratio.md)
  : The Inverse Mills Ratio and Its Derivative

- [`mv_derived.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
  : Standard Deviations and Correlations of a Multivariate Gaussian

- [`mv_derived.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvStudentTDistrib.md)
  : Scale Standard Deviations and Correlations of a Multivariate t

- [`mv_derived.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/mv_derived.multivariate_distrib.md)
  : Matrix Entries as the Default Interpretable Quantities

- [`mv_entry_index()`](https://statmodels7.github.io/distributions7/reference/mv_entry_index.md)
  : Distinct Entries of a Symmetric Matrix, and Their Labels

- [`mv_flat_theta()`](https://statmodels7.github.io/distributions7/reference/mv_flat_theta.md)
  : Require Scalar Parameters

- [`mv_hess_indices()`](https://statmodels7.github.io/distributions7/reference/mv_hess_indices.md)
  : Index Pairs Behind the Hessian Keys of a Multivariate Distribution

- [`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md)
  : The First p Parameters, Read as a Location

- [`mv_location.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md)
  : Mean of a Multivariate Gaussian

- [`mv_location.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/mv_location.MvStudentTDistrib.md)
  : Location of a Multivariate Student t

- [`mv_location.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/mv_location.multivariate_distrib.md)
  : No Location Without a Family That Has One

- [`mv_marginal.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvGaussianDistrib.md)
  : Marginal of a Multivariate Gaussian

- [`mv_marginal.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md)
  : Marginal of a Multivariate Student t

- [`mv_marginal.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/mv_marginal.multivariate_distrib.md)
  : No Marginal Without a Closed Form

- [`mv_moment_start()`](https://statmodels7.github.io/distributions7/reference/mv_moment_start.md)
  : The Moment Estimates a Multivariate Family Starts From

- [`mv_pairs_panels()`](https://statmodels7.github.io/distributions7/reference/mv_pairs_panels.md)
  : Draw the Panel Matrix of a Multivariate Density

- [`mv_prefixed_names()`](https://statmodels7.github.io/distributions7/reference/mv_prefixed_names.md)
  : Prefix a Structure's Free Names with the Matrix They Describe

- [`mv_refuse()`](https://statmodels7.github.io/distributions7/reference/mv_refuse.md)
  : Refuse a Quantity That Has No Multivariate Counterpart

- [`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
  : Standard Deviations and Correlations of a Structured Matrix

- [`mv_sigma.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
  : The Covariance a Multivariate Gaussian Carries

- [`mv_sigma.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvStudentTDistrib.md)
  : The Scale Matrix of a Multivariate Student t

- [`mv_sigma_derivs()`](https://statmodels7.github.io/distributions7/reference/mv_sigma_derivs.md)
  : Derivatives of the Covariance with Respect to Every Parameter

- [`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
  : The Pieces a Multivariate Gaussian Evaluates From

- [`mvg_residuals()`](https://statmodels7.github.io/distributions7/reference/mvg_residuals.md)
  : Residuals and Whitened Residuals

- [`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md)
  : The Pieces a Multivariate t Evaluates From

- [`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
  : The Weight a Multivariate t Gives Each Observation

- [`n_support_points()`](https://statmodels7.github.io/distributions7/reference/n_support_points.md)
  : Number of Points in a Discrete Support

- [`new_check()`](https://statmodels7.github.io/distributions7/reference/new_check.md)
  : Record One Check Result

- [`numDeriv_grad()`](https://statmodels7.github.io/distributions7/reference/numDeriv_grad.md)
  : A Central-Difference Gradient Without a Dependency

- [`observed_deriv()`](https://statmodels7.github.io/distributions7/reference/observed_deriv.md)
  : Observed Derivatives of a Given Order

- [`order_indices()`](https://statmodels7.github.io/distributions7/reference/order_indices.md)
  : Multi-Indices of a Given Order, as Parameter Names

- [`owen_t()`](https://statmodels7.github.io/distributions7/reference/owen_t.md)
  : Owen's T Function

- [`parent_ell()`](https://statmodels7.github.io/distributions7/reference/parent_ell.md)
  : Look Up the Parent's Derivative Components by Block

- [`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md)
  : Probability the Parent Puts on a Single Point

- [`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
  : CDF Gradient When Only Some Parameters Are Location-Scale

- [`plot.continuous_distrib`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
  : Plot Method for Continuous Distributions

- [`plot.discrete_distrib`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md)
  : Plot Method for Discrete Distributions

- [`plot.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
  : Panels of a Multivariate Density

- [`print.FisherScoring`](https://statmodels7.github.io/distributions7/reference/print.FisherScoring.md)
  : Print a Fisher Scoring Specification

- [`print.distrib`](https://statmodels7.github.io/distributions7/reference/print.distrib.md)
  :

  Print Method for `distrib` Objects

- [`print_check_table()`](https://statmodels7.github.io/distributions7/reference/print_check_table.md)
  : Print a Validation Table

- [`safe_check()`](https://statmodels7.github.io/distributions7/reference/safe_check.md)
  : Run a Check, Turning an Error Into a Failure

- [`set_partitions()`](https://statmodels7.github.io/distributions7/reference/set_partitions.md)
  : All Set Partitions of a Finite Index Set

- [`skewness.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.GumbelDistrib.md)
  : Skewness of the Gumbel Distribution

- [`skewness.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.LaplaceDistrib.md)
  : Skewness of the Laplace Distribution

- [`skewness.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.NegBinDistrib.md)
  : Skewness of the Negative Binomial Distribution

- [`skewness.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.PseudoHuberDistrib.md)
  : Skewness of the Pseudo-Huber Distribution

- [`skewness.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormalDistrib.md)
  : Skewness of the Skew Normal Distribution

- [`skewness.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.SkewTDistrib.md)
  : Skewness of the Skew t Distribution

- [`skewness.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/skewness.WeibullDistrib.md)
  : Skewness of the Weibull Distribution

- [`skewness.distrib`](https://statmodels7.github.io/distributions7/reference/skewness.distrib.md)
  : Skewness of a Distribution

- [`skewness.numeric`](https://statmodels7.github.io/distributions7/reference/skewness.numeric.md)
  : Sample Skewness

- [`skewnormal_delta()`](https://statmodels7.github.io/distributions7/reference/skewnormal_delta.md)
  : The Shape a Skew Normal's Moments Depend On

- [`skewt_moment_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_moment_pieces.md)
  : The Quantities a Skew t's Moments Are Built From

- [`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md)
  : The Step a Skew t Differences the Degrees of Freedom With

- [`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
  : The Pieces a Skew t Evaluates From

- [`split_index()`](https://statmodels7.github.io/distributions7/reference/split_index.md)
  : Split a Multi-Index Into Parent and Wrapper Parts

- [`split_mix_theta()`](https://statmodels7.github.io/distributions7/reference/split_mix_theta.md)
  : Split a Wrapper's Parameters From Its Parent's

- [`std_dev.distrib`](https://statmodels7.github.io/distributions7/reference/std_dev.distrib.md)
  : Standard Deviation of a Distribution

- [`std_dev.numeric`](https://statmodels7.github.io/distributions7/reference/std_dev.numeric.md)
  : Sample Standard Deviation

- [`struct_free_or_fit()`](https://statmodels7.github.io/distributions7/reference/struct_free_or_fit.md)
  : Project a Matrix onto What a Structure Can Represent

- [`struct_pair_lookup()`](https://statmodels7.github.io/distributions7/reference/struct_pair_lookup.md)
  : Where Each Pair of Free Values Sits in a Structure's Second
  Derivatives

- [`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
  : Convert Parameter-Scale Derivatives to the Link Scale

- [`trans_deriv_k()`](https://statmodels7.github.io/distributions7/reference/trans_deriv_k.md)
  : Derivatives of a Transformed Distribution

- [`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
  : Second-Order Truncated Moment of the Parent's Derivatives

- [`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
  : Distribution Function of a Truncated Distribution

- [`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
  : The Truncation Constant and Lower Tail

- [`trunc_deriv_k()`](https://statmodels7.github.io/distributions7/reference/trunc_deriv_k.md)
  : Derivatives of a Truncated Distribution

- [`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
  : Expected Hessian of a Truncated Distribution

- [`trunc_gradient()`](https://statmodels7.github.io/distributions7/reference/trunc_gradient.md)
  : Score of a Truncated Distribution

- [`trunc_hess_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_hess_mean.md)
  : Truncated Mean of the Parent's Hessian

- [`trunc_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_hessian.md)
  : Observed Hessian of a Truncated Distribution

- [`trunc_inside()`](https://statmodels7.github.io/distributions7/reference/trunc_inside.md)
  : Which Observations Lie in the Truncated Support

- [`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
  : Derivatives of the Truncation Constant via the Parent's CDF

- [`trunc_pdf()`](https://statmodels7.github.io/distributions7/reference/trunc_pdf.md)
  : Density of a Truncated Distribution

- [`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md)
  : Quantile Function of a Truncated Distribution

- [`trunc_rng()`](https://statmodels7.github.io/distributions7/reference/trunc_rng.md)
  : Random Generation From a Truncated Distribution

- [`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
  : Mean of the Parent's Score Under the Truncated Law

- [`trunc_score_mean_quad()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean_quad.md)
  : Truncated Score Mean by Quadrature

- [`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
  : Truncated Mean of Products of Scores

- [`trunc_y_deriv()`](https://statmodels7.github.io/distributions7/reference/trunc_y_deriv.md)
  : Response Derivative of a Truncated Distribution

- [`variance.GumbelDistrib`](https://statmodels7.github.io/distributions7/reference/variance.GumbelDistrib.md)
  : Variance of the Gumbel Distribution

- [`variance.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md)
  : Variance of the Laplace Distribution

- [`variance.MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/variance.MvGaussianDistrib.md)
  : Variance of a Multivariate Gaussian

- [`variance.MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
  : Covariance of a Multivariate Student t

- [`variance.NegBinDistrib`](https://statmodels7.github.io/distributions7/reference/variance.NegBinDistrib.md)
  : Variance of the Negative Binomial Distribution

- [`variance.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/variance.PseudoHuberDistrib.md)
  : Variance of the Pseudo-Huber Distribution

- [`variance.SkewNormalDistrib`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormalDistrib.md)
  : Variance of the Skew Normal Distribution

- [`variance.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/variance.SkewTDistrib.md)
  : Variance of the Skew t Distribution

- [`variance.WeibullDistrib`](https://statmodels7.github.io/distributions7/reference/variance.WeibullDistrib.md)
  : Variance of the Weibull Distribution

- [`variance.distrib`](https://statmodels7.github.io/distributions7/reference/variance.distrib.md)
  : Variance of a Distribution

- [`variance.numeric`](https://statmodels7.github.io/distributions7/reference/variance.numeric.md)
  : Sample Variance

- [`weibull_gamma_factors()`](https://statmodels7.github.io/distributions7/reference/weibull_gamma_factors.md)
  : Gamma Factors of a Weibull's Moments

- [`weibull_pieces()`](https://statmodels7.github.io/distributions7/reference/weibull_pieces.md)
  : The Pieces a Weibull Evaluates From

- [`za_cont_deriv_k()`](https://statmodels7.github.io/distributions7/reference/za_cont_deriv_k.md)
  : Derivatives of a Zero-Adjusted Continuous Distribution

- [`za_disc_deriv_k()`](https://statmodels7.github.io/distributions7/reference/za_disc_deriv_k.md)
  : Derivatives of a Zero-Adjusted Discrete Distribution

- [`za_y_deriv()`](https://statmodels7.github.io/distributions7/reference/za_y_deriv.md)
  : Response Derivative of a Zero-Adjusted Distribution

- [`zi_deriv_k()`](https://statmodels7.github.io/distributions7/reference/zi_deriv_k.md)
  : Derivatives of a Zero-Inflated Distribution
