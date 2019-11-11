#' @title Sample Expected Eigenvalues from Marcenko-Pastur distribution
#'
#' @description This function sample scaled expected eigenvalues from Marcenko-Pastur distribution with package "RMTsata"(https://cran.r-project.org/web/packages/RMTstat/index.html).
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param params A list of ouput from function CheckDimMatrix.
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
#'   \item{eigenvec:}{ right singular vectors of X matrix with truncation up to dimension rnk.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1,orthogonl = FALSE)
#' MPSamples<-MarcenkoPasturSample(X,rnk=40,times=100)
#' 
#' #equivelantly, if CheckDimMatrix is calcualted
#' params <- CheckDimMatrix(X,rnk=40)
#' MPSamples = MarcenkoPasturSample(params=params,times=100)
#' }
#' should import RMTstat after fix bug
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @export

MarcenkoPasturSample <- 
function(X,                   #data matrix
        params=NULL,          #A list of ouput from function CheckDimMatrix.
        rnk = NA,             #number of singular vectors to estimate
        times=NA,              #split data into X times for parallel computation.
        ...) 
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------

  if(!missing(params)){
    cat("parameters have already been calculated.\n")
  }else{
    if (is.null(X)){
      stop("Invalid input X")
      }else{
        if (missing(rnk)){
        rnk <- min(nrow(X),ncol(X))
        cat('No rnk specified. Calculating full singular value decomposition instead.new rnk = ',rnk,'\n')
      }
        if(missing(times)) times <- 0
        else if (times < 0) stop("times must be positive")
        else if (times > min(nrow(X) - 1, ncol(X) - 1)) stop("times must be strictly less than min(nrow(X), ncol(X))")
  
        params <- CheckDimMatrix(X,rnk=rnk)
        cat("finish checking dimension of X and calculating eigenvalues of X.\n")
      }
  }
  if(missing(times)) times <- 0
  else if (times < 0) stop("times must be positive")
  else if (times > min(params$ndf - 1, params$pdim - 1)) stop("times must be strictly less than min(nrow(X), ncol(X))")
   
  ndf = params$ndf
  pdim = params$pdim
  svr = params$svr
  if(!missing(rnk)&&rnk!=params$rnk) warning("rnk specified does not match rnk calculated in params, use rnk in params instead.\n")
  rnk = params$rnk
  irl = params$irl
  transpose_flag = params$transpose_flag

  var_correct=min(irl$eigen)/MarchenkoPasturPar(ndf,pdim,var=1,svr=svr)$upper#0.00061704 001C
  cat("The corrected variance of MP distribution is ",var_correct,".\n",sep="")
  if(!transpose_flag&&times==0){
    system.time({sim=rmp(pdim, ndf=ndf, pdim=pdim, var=var_correct, svr=ndf/pdim)})
    }else{
      ncores=8
      registerDoParallel(ncores) 
      system.time({
        sim <- foreach(1:times, .combine=c) %dopar% {
          tmp <- rmp(pdim/times, ndf=ndf, pdim=pdim, var=var_correct, svr=ndf/pdim)
        }})
        }
  cat("finish MP sampling.\n")
  MP_irl = tibble(eigen = sim[order(sim,decreasing = T)][1:rnk], dim = 1:rnk)
  
  return(list(ndf=ndf,pdim=pdim,var_correct=var_correct,rnk=rnk,irl=irl,MP_irl=MP_irl,eigenvec = params$eigenvec))
}
