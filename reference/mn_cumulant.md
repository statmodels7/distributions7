# A Joint Cumulant of the Indicators of a Categorical Law

Computes the joint cumulant of the indicators \\1\\Z = t_1\\, \ldots,
1\\Z = t_r\\\\ of a categorical variable with probabilities `prob`,
which is the derivative of the log-sum-exp \\\mathrm{LSE}(\eta)\\ in
\\\eta\_{t_1}, \ldots, \eta\_{t_r}\\. It is the sum over the set
partitions \\\pi\\ of the positions of
\\(-1)^{\lvert\pi\rvert-1}(\lvert\pi\rvert-1)!\prod\_{B\in\pi}\mu_B\\,
where \\\mu_B = p_k\\ if every index of \\B\\ is \\k\\ and zero
otherwise.

## Usage

``` r
mn_cumulant(prob, t)
```

## Arguments

- prob:

  The probability vector.

- t:

  An integer vector of category indices, of length 1 to 4.

## Value

A single number.

## See also

[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md),
the consumer.

## Examples

``` r
p <- c(0.2, 0.3, 0.5)
# The second cumulant is p_k (delta_kl - p_l).
c(distributions7:::mn_cumulant(p, c(1, 1)), 0.2 * 0.8)
#> [1] 0.16 0.16
c(distributions7:::mn_cumulant(p, c(1, 2)), -0.2 * 0.3)
#> [1] -0.06 -0.06
```
