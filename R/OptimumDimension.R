#' @title Rank Estimation in High-Dimensional matrix X.
#'
#' @description Estimate a low rank approximation of a signal-rich subspace in large high-dimensional data.
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param MPSamples A list. Samples from Marc\u{e}nko-Pastur (MP) distribution and eigenvalues of X. Output from MarcenkoPasturSample function.
#' @param rnk number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @param times split data into X times for parallel computation.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ number of degrees of freedom of X and the white Wishart matrix.}
#'   \item{pdim:}{ number of dimensions of X and the white Wishart matrix.}
#'   \item{var_correct:}{ population variance for Marcenko-Pastur distribution.}
#'   \item{irl:}{ a data frame of scaled eigenvalues and corresponding dimensions.}
#'   \item{MP_irl:}{ a data frame of samped expected eigenvalues from Marcenko-Pastur and corresponding dimensions.}
#'   \item{bcp_irl:}{ probability of change in mean and posterior means of eigenvalue difference between $X$ and $N$.}
#'   \item{changePoint:}{ estimated changepoint position by "cpm" package.}
#' }
#' @section Details:
#' We estimate a low rank approximation of a signal-rich subspace in large high-dimensional data by decomposing matrix
#' into a signal-plus-noise space and approximate the signal-rich subspace with a rank K approximation
#' \eqn{\hat{X}=\sum_{k=1}^{K}d_ku_k{v_k}^T}. To estimate rank K, we propose a simple procedure assuming that matrix X is composed
#' of a low-rank signal matrix S and an average general noise random matrix \eqn{\bar{N}}. It has been shown that
#' the average eigenvalues of random matrices N follows a universal Marc\u{e}nko-Pastur (MP) distribution.
#' We hypothesize that the deviation of eigenvalues of X from the MP distribution indicates the intrinsic dimension of signal-rich subspace.
#' @examples
#' \donttest{
#' results<-OptimumDimension(X,rnk=10,times=1000)
#' str(results)
#' ScreePlot(results$MarcenkoPasturSample,annotation=0)
#' modified_legacyplot(results$Changepoint$bcp.irl,annotation=10)
#' }
#' should import RMTstat after fix bug
#' @seealso \code{\link[RMTstat]} for details of Marcenko-Pastur distribution.
#' @importFrom bcp bcp
#' @importFrom cpm detectChangePoint
#' @importFrom  tibble tibble
#' @export

OptimumDimension <- 
  function(X,                   #data matrix
        MPSamples=NULL,          #A list of ouput from function MarcenkoPasturSample.
        rnk = NA,             #number of singular vectors to estimate
        times=NA,              #split data into X times for parallel computation.
        ...) {


  if (is.null(X))  stop("Invalid input X")
  if(missing(times)) {times <- 0}
  else if (times < 0)  stop("times must be positive")
  else if (times > min(nrow(X) - 1, ncol(X) - 1)) stop("times must be strictly less than min(nrow(X), ncol(X))")
  if (!missing(MPSamples)) cat("MP samples have already been calculated.\n")
  if (missing(rnk)) {rnk <- min(nrow(X),ncol(X));cat('No rnk specified. Calculating full singular value decomposition instead.\n');cat('rnk missing, new rnk = ',rnk,'\n')}
  if (missing(MPSamples)) {MPSamples <- MarcenkoPasturSample(X,rnk=rnk,times=times);cat("finish calculating MP samples.\n")}

  ndf = MPSamples$ndf
  pdim = MPSamples$pdim
  rnk = MPSamples$rnk
  var_correct = MPSamples$var_correct
  irl = MPSamples$irl
  MP_irl = MPSamples$MP_irl

  sigma_a = irl$eigen
  sigma_MP = MP_irl$eigen

  #Bayesian Change Point
  bcp.irl = bcp(as.vector(sigma_a-sigma_MP), p0 = 0.1)
  #Bayesian Posterior Prob Change Point
  bcp_post = bcp(as.vector(c(bcp.irl$posterior.prob[-rnk],0)),p0=0.1)
  #Bayesian Diff Change Point
  bcp_diff = bcp(as.vector(-diff(sigma_a-sigma_MP)),p0=0.001)
  #Single Change Point
  changePoint = detectChangePoint(sigma_a-sigma_MP,cpmType="Exponential")$changePoint
  
  return(list(MarcenkoPasturSample=list(ndf=ndf,pdim=pdim,rnk=rnk,var_correct=var_correct,irl=irl,MP_irl=MP_irl),
              Changepoint=list(bcp.irl=bcp.irl,bcp_post=bcp_post,bcp_diff=bcp_diff,changePoint=changePoint)))
}