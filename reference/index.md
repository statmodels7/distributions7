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
- [`gamma_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
  : Gamma Distribution Object (Mean-Variance Parameterization)
- [`invgauss_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
  : Inverse-Gaussian Distribution Object (Mean-Dispersion
  Parameterization)
- [`lognormal_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
  : Lognormal Distribution Object (Log-Scale Parameterization)
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
- [`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
  : Derivatives on the Link (Real) Scale
- [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  : Strategies for Expected Derivatives

## Estimation

- [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  : Maximum-Likelihood Estimation
- [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  : S7 Class for Maximum-Likelihood Fits
- [`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
  : Numerically Validate a Distribution

## Moments

- [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  : Calculate the Expected Value of a Function
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
point mass at zero, a restriction to an interval, or the law of a
transformed variable.

- [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  : Zero-Inflated Distribution Object (Discrete)
- [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
  : Zero-Adjusted Distribution Object
- [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
  : Truncated Distribution Object
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

## Classes

The S7 classes; each page lists the methods that dispatch on it.

- [`distrib()`](https://statmodels7.github.io/distributions7/reference/distrib.md)
  : S7 Class for Probability Distributions
- [`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
  : S7 Class for Continuous Distributions
- [`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
  : S7 Class for Discrete Distributions

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
- [`distrib_cdf.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GaussianDistrib.md)
  : Gaussian Cumulative Distribution Function
- [`distrib_deriv3.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GaussianDistrib.md)
  : Gaussian Analytical Third-Order Derivatives
- [`distrib_deriv4.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GaussianDistrib.md)
  : Gaussian Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GaussianDistrib.md)
  : Gaussian Analytical Expected Hessian
- [`distrib_grad_y.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GaussianDistrib.md)
  : Gaussian Response Derivatives
- [`distrib_gradient.GaussianDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GaussianDistrib.md)
  : Gaussian Analytical Gradient
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
- [`distrib_grad_y.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.CauchyDistrib.md)
  : Cauchy Response Derivatives
- [`distrib_gradient.CauchyDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
  : Cauchy Analytical Gradient
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
- [`distrib_grad_y.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LogisticDistrib.md)
  : Logistic Response Derivatives
- [`distrib_gradient.LogisticDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
  : Logistic Analytical Gradient
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
- [`distrib_cdf.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentTDistrib.md)
  : Student's t Cumulative Distribution Function
- [`distrib_deriv3.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentTDistrib.md)
  : Student's t Analytical Third-Order Derivatives
- [`distrib_deriv4.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentTDistrib.md)
  : Student's t Analytical Fourth-Order Derivatives
- [`distrib_expected_hessian.StudentTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentTDistrib.md)
  : Student's t Analytical Expected Hessian
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
- [`distrib_expected_hessian.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
  : Laplace Analytical Expected Hessian (Fisher Information)
- [`distrib_grad_y.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md)
  : Laplace Response Derivatives
- [`distrib_gradient.LaplaceDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md)
  : Laplace Analytical Gradient
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

## Wrapper classes

The methods of the zero-inflated, zero-adjusted, truncated and
transformed wrappers.

- [`ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md)
  : S7 Class for Zero-Inflated Distributions
- [`distrib_cdf.ZeroInflatedDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroInflatedDistrib.md)
  : Zero-Inflated Cumulative Distribution Function
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
- [`invgauss_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
  : Inverse-Gaussian Distribution Object (Mean-Dispersion
  Parameterization)
- [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
  : Laplace Distribution Object
- [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
  : Logistic Distribution Object
- [`lognormal_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
  : Lognormal Distribution Object (Log-Scale Parameterization)
- [`negbin_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
  : Negative Binomial Distribution Object (NB2)
- [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
  : Poisson Distribution Object
- [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
  : Pseudo-Huber Distribution Object (Location-Scale Parameterization)
- [`student_t_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t_distrib.md)
  : Student's t Distribution Object (Location-Scale Parameterization)
