library(gridExtra)

extract_legend <- function(p) {
  gt <- ggplot_gtable(ggplot_build(p))
  legidx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  legend <- gt$grobs[[legidx]]
  return(legend)
}