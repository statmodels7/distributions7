# The three ready-made strategies that fit_distrib(method = ) took as strings,
# written as the objects that replace them.
named_methods <- function() {
  list(
    fisher = fisher_scoring(),
    newton = optimizers7::newton(criterion = optimizers7::crit_grad()),
    bfgs = optimizers7::bfgs(criterion = optimizers7::crit_grad())
  )
}
