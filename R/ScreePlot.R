#' @title Scree Plot of scaled eigenvalues of X and sampled Expected Eigenvalues from Marcenko-Pastur distribution
#'
#' @description This function plot scree plot of scaled eigenvalues of X and sampled scaled expected eigenvalues from Marcenko-Pastur distribution.
#' @param X A list. Samples from Marc\u{e}nko-Pastur (MP) distribution and eigenvalues of X. Output from MarcenkoPasturSample function.
#' @examples
#' \donttest{
#' ScreePlot(MPSamples,Changepoint=2,annotation=0)
#' }
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom ggrepel geom_text_repel
#' @export

ScreePlot <- function(MPSamples,Changepoint=NULL,annotation=NULL) {

  rnk= nrow(MPSamples$irl)
  ndf = MPSamples$ndf
  pdim = MPSamples$pdim
  var_correct = MPSamples$var_correct
  transpose_flag = MPSamples$transpose_flag
  irl = MPSamples$irl
  MP_irl = MPSamples$MP_irl
  if (missing(annotation)) {annotation <- rnk;cat("Anotating for all point. Set annotation = 0 to stop annotation.\n")}
  if (!is.null(annotation)&annotation > rnk) stop("Annotation number must be strictly less or equal to than rnk")
  if(annotation==0) {mark=rep("",rnk)}else{
  mark=c(1:annotation,rep("",rnk-annotation))}

  noise_scree = ggplot() + geom_line(aes(x = dim, y = eigen),irl,colour="black") + geom_point(aes(x = dim, y = eigen),irl,color = "red") +
    geom_line(aes(x = dim, y = eigen),MP_irl,colour="black") + geom_point(aes(x = dim, y = eigen),MP_irl,color = "blue") +
    theme_minimal() + xlab("Dimension") + ylab("Eigenvalue Scaled")+theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Scree Plot","\n","N = ",ifelse(transpose_flag,pdim,ndf),", P = ",ifelse(transpose_flag,ndf,pdim),", Var = ",round(var_correct,2))) +
    geom_text_repel(aes(x= irl$dim, y = irl$eigen, label=mark),colour="black",size=5) 

  if(!missing(Changepoint)) noise_scree = noise_scree + ggtitle(paste0("Scree Plot","\n","N = ",ifelse(transpose_flag,pdim,ndf),", P = ",ifelse(transpose_flag,ndf,pdim),", Var = ",round(var_correct,2),", ChnagePoint est = ",Changepoint))
  return(noise_scree)
}


