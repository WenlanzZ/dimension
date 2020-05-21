#' Modified K means clustering
#'
#' k-means clustering modified to cluster probabilities of change points
#'  to 2 clusters and return total within-cluster variation. Adapted from 
#' https://github.com/TunChiehHsu/Kmeans.git.
#'
#' @param data     an n by p numeric matrix; the data
#'                matrix of predictors
#' @param k     integer giving the number of neighbors to
#'              include. k = 2.
#'
#' @return      a data frame of total within-cluster variation.

km <- function(data) {
  prev_center <- NULL
  update_center <- NULL
  it <- 1
  k <- 2
  within_var <- matrix(0, nrow = nrow(data), ncol = 3)
  while (it < nrow(data) + 1) {
    if (it == 1) {
      center <- list(data[1, ], data[nrow(data), ])
    }
    #record distance
    gd <- NULL
    for (i in 1:k) {
      # calculate distance based on different center
      temp_center <- center[[i]]
      #distance function
      dist_m <- function(a) {
        sum((a - temp_center)^2)
      }
      #calculate distance
      group <- apply(data, 1, function(x) dist_m(x))
      #record
      gd <- cbind(gd, group)
    }
    gd <- data.frame(gd)
    colnames(gd) <- paste("group", c(1:k))
    gd$label <- c(rep(1, it), rep(2, nrow(data) - it))

    # calculate total within-cluster variation
    g1_within <- sum(gd[which(gd$label == 1), 1])
    g2_within <- sum(gd[which(gd$label == 2), 2])
    tol_within <- sum(g1_within, g2_within)
    within_var[it, ] <- c(tol_within, g1_within, g2_within)
    prev_center <- center
    center <- list()
    # calculate new center
    for (i in 1:k) {
      sub_data <- data[which(gd$label == i), ]
      update_center <- colMeans(sub_data)
      center[[i]] <- update_center
      sub_data <- NULL
    }
    update_center <- center
    it <- it + 1
  }
  return(within_var)
}
