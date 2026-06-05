# ConvRt :: plot.R
#
# Convenience ggplot2 helper for a fitted R_t curve + confidence band.

# Plot an R_t estimate with its confidence band.
#
#   x      : a fit_convrt_* result (uses its $rt_df), or an rt_df tibble directly
#            (columns date, Rt_mean, Rt_lo, Rt_hi).
#   truth  : optional tibble(date, Rt_true) to overlay (e.g. in simulations).
#   title  : optional plot title.
#
# Returns a ggplot object (so you can add layers / theming).  Requires ggplot2.
plot_convrt <- function(x, truth = NULL, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_convrt() requires the ggplot2 package. install.packages(\"ggplot2\")",
         call. = FALSE)
  rt_df <- if (is.data.frame(x)) x else x$rt_df
  need  <- c("date", "Rt_mean", "Rt_lo", "Rt_hi")
  if (is.null(rt_df) || !all(need %in% names(rt_df)))
    stop("plot_convrt(): expected a fit_convrt_* result or an rt_df with columns ",
         paste(need, collapse = ", "), ".")
  rt_df$date <- as.Date(rt_df$date)

  p <- ggplot2::ggplot(rt_df, ggplot2::aes(date)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Rt_lo, ymax = Rt_hi),
                         fill = "#0072B2", alpha = 0.22) +
    ggplot2::geom_line(ggplot2::aes(y = Rt_mean), colour = "#0072B2", linewidth = 1.0)
  if (!is.null(truth)) {
    truth$date <- as.Date(truth$date)
    p <- p + ggplot2::geom_line(data = truth,
                                ggplot2::aes(date, Rt_true),
                                colour = "grey10", linewidth = 1.1)
  }
  p +
    ggplot2::labs(title = title, x = NULL, y = expression(R[t])) +
    ggplot2::theme_bw(base_size = 14)
}
