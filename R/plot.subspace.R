#' @title Scree plot of scaled eigenvalues of X and random noise matrix N
#'
#' @description This is a generic function for supspace class to plot scree plot of scaled eigenvalues of X and sampled scaled expected eigenvalues from Marcenko-Pastur distribution.
#' @param data A subsapce class.
#' @param Changepoint A number. Estimated changepoint in OptimumDimension function.
#' @param annotation A number. Choose to label points up to annotation number. Set to 0 with no annotation.
#' @examples
#' \donttest{
#' plot(Subspace, Changepoint = 0, annotation = 15)
#' }
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom ggrepel geom_text_repel
#' @export

plot.subspace <- function(obj,                              # A subspace class
                          Changepoint = NULL,               # Estimated changepoint in OptimumDimension function
                          annotation = NULL,                # Choose to label points up to annotation number
                          verbose = TRUE)
{
# ---------------------------------------------------------------------------------------------------------
# Basic parameter set up
# ---------------------------------------------------------------------------------------------------------
  ndf             <- obj$ndf
  pdim            <- obj$pdim
  rank            <- obj$rank
  rnk             <- max(rank)
  var_correct     <- obj$var_correct
  transpose_flag  <- obj$transpose_flag
  irl             <- obj$irl
  MP_irl          <- obj$MP_irl
# ---------------------------------------------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------------------------------------------
  if (missing(annotation)) {
    annotation <- rnk
    if (verbose) {
      cat("Anotating for all points. Set annotation = 0 to stop annotation.\n")
    }
  }
  if (!is.null(annotation) & annotation > rnk) 
  {
    annotation <- rnk
    warning("Annotation number must be strictly less or equal to than maximum rank.\n")
  }
  if (annotation == 0) {
    mark <- rep("", length(rank))
  } else {
    mark <- ifelse(annotation >= rank, rank, "")
  }
# ---------------------------------------------------------------------------------------------------------
# Scree plot for both eigenvalues of X and simulated eigenvalues of noise
# ---------------------------------------------------------------------------------------------------------
  scree <- ggplot() + 
            geom_line(aes(x = dim, y = eigen), irl, colour = "black") + 
            geom_point(aes(x = dim, y = eigen), irl, color = "red") +
            theme_minimal() + 
            xlab("Dimension") + 
            ylab("Eigenvalue Scaled") + 
            theme(plot.title = element_text(hjust = 0.5)) + 
            geom_text_repel(aes(x= irl$dim, y = irl$eigen, label=mark), colour="black", size=5) +
            ggtitle(
              paste0("Scree Plot\n",
                     "N = ", 
                     ifelse(transpose_flag,pdim,ndf), 
                     ", P = ", 
                     ifelse(transpose_flag,ndf,pdim), 
                     ", Var = ", 
                     round(var_correct,2))
            )

  if (!is.null(MP_irl)) {
    scree <- scree +           
              geom_line(aes(x = dim, y = eigen), MP_irl, colour = "black") + 
              geom_point(aes(x = dim, y = eigen), MP_irl, color = "blue")
  }
          
  if (!missing(Changepoint)) {
    scree <- scree + 
              ggtitle(
                paste0("Scree Plot\n", 
                       "N = ", 
                       ifelse(transpose_flag,pdim,ndf), 
                       ", P = ", 
                       ifelse(transpose_flag,ndf,pdim), 
                       ", Var = ", 
                       round(var_correct,2), 
                       ", ChnagePoint est = ", 
                       Changepoint)
              )
  }

  return(scree)
}
