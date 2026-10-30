#' The autoplot() Generic of ggplot2
#'
#' @description
#' The generic `autoplot()` belongs to \pkg{ggplot2}, which this package
#' suggests and does not import. The methods written on it here are
#' registered by `S7::methods_register()` when \pkg{ggplot2} is loaded.
#'
#' @return An S7 external generic, the description of `ggplot2::autoplot()`
#'   that `S7::method()` registers methods on.
#'
#' @keywords internal
autoplot <- S7::new_external_generic("ggplot2", "autoplot", "object")


#' The Packages That a Graph Requires
#'
#' @description
#' Signals an error that names the missing packages when \pkg{ggplot2}, or
#' \pkg{patchwork} with `patchwork = TRUE`, is not installed. The graphs of
#' the package are drawn with \pkg{ggplot2}, which is suggested and not
#' imported.
#'
#' @param patchwork Whether the graph also requires \pkg{patchwork}.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
need_ggplot <- function(patchwork = FALSE) {
  pkgs <- c("ggplot2", if (patchwork) "patchwork")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(sprintf("The graph is drawn with %s; install %s.",
                 paste(miss, collapse = " and "),
                 if (length(miss) == 1L) "it" else "them"), call. = FALSE)
  }
  invisible(NULL)
}


#' The Arguments That a Graph Accepts
#'
#' @description
#' Signals an error when `...` is not empty, and when `legend` is not `TRUE`
#' or `FALSE`. The graphs take no graphical arguments: the returned ggplot2
#' object is modified instead.
#'
#' @param legend The value of the argument `legend`.
#' @param ... The `...` of the calling method.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
check_plot_args <- function(legend, ...) {
  if (...length()) {
    stop("The graph takes no further arguments; modify the returned object ",
         "instead.", call. = FALSE)
  }
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("'legend' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}


#' The Placement of the Legends of a Graph
#'
#' @description
#' Places the legends below the panel, one above the other, or removes them
#' when `legend` is `FALSE`.
#'
#' @param legend Whether the graph draws its legends.
#'
#' @return A `ggplot2::theme()` object.
#'
#' @keywords internal
plot_key_theme <- function(legend = TRUE) {
  if (!legend) return(ggplot2::theme(legend.position = "none"))
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                 legend.spacing.y = ggplot2::unit(2, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))
}


#' Whole-Number Breaks of an Axis
#'
#' @description
#' Returns the breaks that `pretty()` chooses over `limits`, keeping the
#' whole numbers, for the horizontal axis of a graph over the support of a
#' discrete family.
#'
#' @param limits The limits of the axis, a numeric vector of length 2.
#'
#' @return A numeric vector of whole numbers.
#'
#' @keywords internal
integer_breaks <- function(limits) {
  b <- pretty(limits)
  b[b == round(b)]
}
