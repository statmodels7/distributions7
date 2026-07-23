# Default Numerical RNG for Discrete Distributions

Fallback method: discrete distributions that do not implement a native
RNG generate draws by inverse transform sampling on the cumulative pmf
table.

No rejection scheme is needed here, and none would help. Inverting the
cumulative mass function is exact, because the distribution function of
a lattice-valued variable is a step function: there is nothing to solve.
The table is built once per distinct parameter value and the whole
sample is located in it with a single binary search, which costs a
fraction of a microsecond per draw.

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- n:

  Number of observations to generate.

- theta:

  A named list of parameters.
