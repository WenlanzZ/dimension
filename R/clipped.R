#' @title Eigenvalue clipping Procedure
#'
#' @description This function clips the scaled eigenvalues of X in order to provide a cleaned estimator E_clipped of the underlying correlation matrix. Proceeds by keeping the [N * alpha] top eigenvalues and shrinking the remaining ones by a trace-preserving constant (i.e. Tr(E_clipped) = Tr(E)). This function is adpated from "Python for Random Matrix Theory" credit to J.-P. Bouchaud and M. Potters.
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param params A list of ouput from function CheckDimMatrix.
#' @param rnk number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @param method the method to be used; method = "threshold" returns (1-alpha)*rnk proportion of eigenvalues above threshold; 
#'               method = "hard" returns all the empirical eigenvalues greater than the upper limit of the support to the Marcenko-Pastur spectrum;
#'               method = "identity" returns eigenvalued specified in location vector.
#' @param alpha determining the fraction to keep of the top eigenvalues in threshold method. 
#' @param location indicate the location of eigenvalues to keep in identity method.
#' @param zeroout zero out eigenvalues when clipping. default is to set it to average.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{xi_clipped:}{ cleaned estimator of eigenvalues through a simple eigenvalue clipping procedure (cf. reference below).}
#'   \item{X_clipped:}{ cleaned estimator of the true sigmal matrix underlying a noisy high dimensional matrix X.}
#'   \item{E_clipped:}{ cleaned estimator of the true correlation matrix underlying a noisy, in-sample estimate E (empirical correaltion matrix estimated from X (cf. reference below)).}
#'   \item{v:}{ right singular vectors of X matrix with truncation up to dimension rnk.}
#'   \item{u:}{ left singular vectors of X matrix with truncation up to dimension rnk.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1,orthogonl = FALSE)
#' x_clp<-clipped(x,rnk=20,method="threshold",alpha=0.9,zeroout=TRUE)
#' x_clp<-clipped(x,rnk=20,method="hard",zeroout=FALSE)
#' x_clp<-clipped(x,rnk=20,method="identity",location=c(1:15),zeroout=FALSE)
#' 
#' #equivelantly, if CheckDimMatrix is calcualted
#' params <- CheckDimMatrix(X,rnk=40)
#' x_clp<-clipped(params=params,method="threshold",alpha=0.9,zeroout=TRUE)
#' }
#' Reference
#' ---------
#' "Financial Applications of Random Matrix Theory: a short review",
#' J.-P. Bouchaud and M. Potters
#' arXiv: 0910.1205 [q-fin.ST]
#' should import RMTstat after fix bug
#' @export

clipped <- function(X,                                           # data matrix
                    params=NULL,                                 # a list of ouput from function CheckDimMatrix.
                    rnk=NA,                                      # number of singular vectors to estimate
                    method=c("threshold","hard","identity"),     # choose method from c("threshold","hard","identity")
                    alpha=NA,                                    # determining the fraction to keep of the top eigenvalues of an empirical correlation matrix. 
                    location = NA,                               # indicate the location of eigenvalues to keep.
                    zeroout = FALSE,                             # zero out eigenvalues when clipping. default is to set it to average.
                    ...)                                         # optional arguments
{

# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------

  if(!missing(params)){
    cat("Parameters have already been calculated.\n")
  }else{
    if (is.null(X)){
      stop("Invalid input X")
      }else{
        if (missing(rnk)){
        rnk <- min(nrow(X),ncol(X))
        cat('No rnk specified. Calculating full singular value decomposition instead.new rnk = ',rnk,'\n')
      }
        params <- CheckDimMatrix(X,rnk=rnk)
        cat("Finish checking dimension of X and calculating eigenvalues of X.\n")
      }
  }

    ndf = params$ndf
    pdim = params$pdim
    svr = ndf/pdim
    if(!missing(rnk)&&rnk!=params$rnk) warning("Rnk specified does not match rnk calculated in params, use rnk in params instead.\n")
    rnk = params$rnk
    irl = params$irl
    v = params$v
    u = params$u
    transpose_flag = params$transpose_flag

    switch(method,
        hard={
          lambda_min = MarchenkoPasturPar(ndf,pdim,var=1,svr=svr)$lower
          cat("lambda_min = ",lambda_min,"\n")
          lambda_max = MarchenkoPasturPar(ndf,pdim,var=1,svr=svr)$upper
          cat("lambda_max = ",lambda_max,"\n")
          xi_clipped = ifelse(irl$eigen>=lambda_max,irl$eigen,NA)
        },
        threshold={
          if(missing(alpha)) stop("Alpha must be specified")
          else if (alpha < 0) stop("Alpha must be positive")
          else if (alpha > 1) stop("Alpha must be less or equal to 1")
          xi_clipped = rep(NA,rnk)
          threshold = ceiling(alpha * rnk)
          cat("threshold = ",threshold,"\n")
          if (threshold > 0) xi_clipped[1:threshold] = irl$eigen[1:threshold]
          else(stop("No eigenvalue preserved"))
        },
        identity={
          if(!is.numeric(location)) stop("Invalid location input")
          if(missing(location)) stop("Location must be specified")
          if(max(location)>rnk) stop ("Location must be smaller than rnk")
          else if(min(location)<=0) stop ("Location out of bounds")
          xi_clipped = rep(NA,rnk)
          xi_clipped[location] = irl$eigen[location]
          },
        stop("Invalid method input")
    )
    gamma = ifelse(zeroout,0,(sum(irl$eigen) - sum(na.omit(xi_clipped)))/sum(is.na(xi_clipped))) 
    cat("gamma = ",gamma,"\n")
    xi_clipped = ifelse(is.na(xi_clipped),gamma,xi_clipped)
    #calculate estimated X
    X_clipped = u%*%diag(xi_clipped*pdim)%*%t(v)
    #calculate empirical covariance
    V_clipped = v%*%diag(xi_clipped*pdim)%*%t(v)/ (ndf - 1L)
    ## symmetric rescaling to correlation matrix
    E_clipped<-V_clipped / tcrossprod(diag(V_clipped) ^ 0.5)
    # ---------------------------------------------------------------------   
    # E_clipped = np.zeros((N, N), dtype=float)
    # for xi, eigvec in zip(xi_clipped, eigvecs):
    #     eigvec = eigvec.reshape(-1, 1)
    #     E_clipped += xi * eigvec.dot(eigvec.T)
        
    # tmp = 1./np.sqrt(np.diag(E_clipped))
    # E_clipped *= tmp
    # E_clipped *= tmp.reshape(-1, 1)
    # ---------------------------------------------------------------------  

    return(list(xi_clipped=xi_clipped,X_clipped=X_clipped,E_clipped=E_clipped,v=v,u=u))
}


