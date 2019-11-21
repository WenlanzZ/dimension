#' @title Return parameter setting for MarcenkoPasturSample function.
#'
#' @param X A numeric real- or real-valued sparse matrix with n number of samples and p number of features. If p>n, a warning message is generated and the transpose of X is used.
#' @param rnk number of right singular vectors to estimate. rnk must be smaller or equal to max(nrow(X),ncol(X)).
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{ndf:}{ number of degrees of freedom of X and the white Wishart matrix.}
#'   \item{pdim:}{ number of dimensions of X and the white Wishart matrix.}
#'   \item{svr:}{ samples to variance ratio; equals to ndf divided by pdim.}
#'   \item{rnk:}{ number of right singular vectors to estimate.}
#'   \item{transpose_flag:}{ whether the matrix X is transposed.}
#'   \item{irl:}{ a data frame of scaled eigenvalues and corresponding dimensions.}
#'   \item{v:}{ right singular vectors of X matrix with truncation up to dimension rnk.}
#'   \item{u:}{ left singular vectors of X matrix with truncation up to dimension rnk.}
#' }
#' @examples
#' \donttest{
#' X <- Xsim(n=1000,p=500,ncc=10,var=2,fact = 1,orthogonl = FALSE)
#' params <- CheckDimMatrix(X,rnk=40)
#' }
#' #should import RMTstat after fix bug
#' @seealso Random Matrix Theory pacakge credit to Gregory Giecold and Lionel Ouaknin
#' @importFrom  tibble tibble
#' @importFrom  irlba irlba
#' @export

CheckDimMatrix <- 
function(X,                 #data matrix
         rnk=NA,             #number of singular vectors to estimate
         ...)               # optional arguments
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------
    if (is.null(X))  stop("Invalid input X")
    ndf = nrow(X)
    pdim = ncol(X)
    if (missing(rnk)) {rnk <- min(ndf,pdim);cat('No rnk specified. Calculating full singular value decomposition instead.\n')}
    if(rnk <= 0) stop("Rnk must be positive")#cat('rnk missing',!missing(rnk),'rnk = ',rnk,'rnk <=0',rnk <= 0,'and',!missing(rnk)&&rnk <= 0,'\n')
    else if(rnk >min(ndf, pdim)) stop("Rnk out of bounds")#cat('rnk missing',!missing(rnk),'rnk = ',rnk,'rnk <=0',rnk <= 0,'and',!missing(rnk)&&rnk <= 0,'\n')
    
    transpose_flag=FALSE
    if(nrow(X)<ncol(X)){
        warning('The number of samples of X is smaller than the number of features of X. A transpose of X is used instead.\n')
        X = t(X)
        transpose_flag=TRUE
        ndf = nrow(X)
        pdim = ncol(X)
    }
    svr = ndf/pdim
    Xstd<-sweep(X, 2L, colMeans(X))
    if(rnk>pdim/2){tmp = svd(Xstd)}else{tmp = irlba(Xstd, nv = rnk)}
    irl = tibble(eigen = tmp$d[1:rnk]^2/(pdim), dim = 1:rnk)
    v = tmp$v[,1:rnk]
    u = tmp$u[,1:rnk]
    
    return(list(ndf=ndf,pdim=pdim,svr=svr,rnk=rnk,transpose_flag=transpose_flag,irl=irl,v=v,u=u))
}
