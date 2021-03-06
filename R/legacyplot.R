#' @title Summary plots of bcp() output
#'
#' @description This function produces summary plots of bcp() output.
#'  It is adapted from the default legacyplot function in pacakge bcp.
#' @param x A list. The result of a call to bcp().
#' @param annotation Annotate points up to total components calculated.
#'  No annotation when annotation = 0.
#' @param medianfilter A logical value. Compute running medians to smooth
#'  scatter plot.
#' @examples
#' library(bcp)
#' legacyplot(bcp(as.vector(c(rep(10, 10), 9.5, rep(0, 10))),
#' p0 = 0.1))
#' @import ggplot2
#' @import gridExtra
#' @importFrom stats runmed
#' @export

legacyplot <- function(x, annotation   = NULL, medianfilter = FALSE) {
  if (class(x) == "dimension") {
    x <- x$bcp_irl
  }
  x$posterior.prob[nrow(x$data)] <- 0
  if (missing(annotation)) {
    annotation <- nrow(x$data)
    message("Anotating for all point. Set annotation = 0 to stop annotation.\n")
  }
  if (!is.numeric(annotation)) {
    stop("Anotation must be numbers.\n")
  }
  if (!is.null(annotation) & max(annotation) > nrow(x$data)) {
    stop("Annotation number must be strictly less or equal to than
         rows of x.\n")
  }
  if (length(annotation) > 1) {
    mark <- rep("", nrow(x$data))
    mark[annotation] <- annotation
  } else {
    mark <- c(1:annotation, rep("", nrow(x$data) - annotation))
    if (annotation == 0) {
      mark <- rep("", nrow(x$data))
    }
  }

  p1 <- ggplot() +
          theme_minimal() +
          xlab("") +
          ylab("Posterior Mean") +
          theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 0.5)) +
          ggtitle(paste0("Posterior Means and Probabilities of a Change\n")) +
          #scale_x_continuous(breaks = seq(0, nrow(x$data), 1)) +
          geom_point(aes(x = x$data[, 1], y = x$data[, 2]),
                     colour = "red") +
          geom_line(aes(x = x$data[, 1], y = x$posterior.mean),
                    color = "black",
                    cex = 0.5) +
          geom_text_repel(aes(x = x$data[, 1], y = x$posterior.mean,
                              label = mark),
                          colour = "black",
                          size = 5)
  p2 <- ggplot() +
          theme_minimal() +
          xlab("Dimension") +
          ylab("Posterior Probability") +
          ylim(0, 1) +
          theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 0.5)) +
          #scale_x_continuous(breaks = seq(0, nrow(x$data), 1)) +
          geom_point(aes(x = x$data[, 1],
                         y = switch(2 - medianfilter,
                                    runmed(x$posterior.prob, 3),
                                    x$posterior.prob)),
                     colour = "blue") +
          geom_line(aes(x = x$data[, 1],
                        y = switch(2 - medianfilter,
                                   runmed(x$posterior.prob, 3),
                                   x$posterior.prob)),
                    color = "black",
                    cex = 0.5) +
          geom_text_repel(aes(x = x$data[, 1],
                              y = switch(2 - medianfilter,
                                         runmed(x$posterior.prob, 3),
                                         x$posterior.prob),
                              label = mark),
                          colour = "black",
                          size = 5)

  #merge all three plots within one grid (and visualize this)
  grid.arrange(p1, p2, ncol = 1)
}
