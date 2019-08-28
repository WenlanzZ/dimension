#' @title Fit the linear exponential model using GenSA
#'
#' @description Fit the eigenvalues of X matrix with linear exponential model.
#' @param X a numeric matrix.
#' @param sacle a logical value.
#' @param d_start the dimention start to fit the linear exponential model.
#' @param rnk number of right singular vectors to estimate.
#' @param lower Vector with length of par. Lower bounds for components.
#' @param upper Vector with length of par. Upper bounds for components.
#' @param plot a logical value. Ouput plots if true.
#' @param times Integer. Total simulation times for noise matrix.
#' @param seed the random number seed.
#' @examples
#' \donttest{
#' # Simulate a matrix X with first two columns correlated.
#' X <- Xsim(n=100,p=90,ncc=2,var=10)
#'
#'
#' res <- main_linear(X,scale=TRUE,d_start=2,rnk=20)
#' }
#' @import foreach
#' @importFrom crayon red
#' @importFrom tibble tibble
#' @importFrom irlba irlba
#' @importFrom data.table data.table
#' @importFrom GenSA GenSA
#' @importFrom stats weighted.mean
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @importFrom ggrepel geom_text_repel
#' @export

#scale = TRUE;d_start = 10;rnk = 20;lower = c(-10,-1000,-100);upper = c(100,1000,100);plot = TRUE;times = 10000;seed = 1000;
main_linear <- function(X,
                        scale = TRUE,
                        d_start = 10,#########------------------------------------------------------->>>>>fix later
                        rnk = 20,#########------------------------------------------------------->>>>>fix later
                        lower = c(-10,-1000,-100),
                        upper = c(100,1000,100),
                        plot = TRUE,
                        times = 10000,
                        seed = 1000) {
  set.seed(seed)
  if (is.null(getDoParName())) {
    registerDoSEQ()
  }
  if (d_start >= p) {
    stop(red(
      "d_start must be less than dimension of X."
    ))
  }

  ### Calculate irlba for matrix X ###
  if(scale){
    E<-scale(X)
  }else{
    E<-X
  }
  n <- nrow(E)
  p <- ncol(E)
  rnk<-round(p/3)#########------------------------------------------------------->>>>>justification for p/3 : laddle paper

  irl <- tibble(sigma = irlba(E, nv = rnk)$d, dim = 1:(rnk))
  names(irl)[1] <-"sigma"

  sq_noise<-foreach(i = 1:times) %dopar% {
    xn<-matrix(rnorm(n*p),n,p)
    xnprime<-xn/sqrt(n)

    irl_noise <- data.table(if(rnk<(p/3)){sigma = irlba(xn, nv = rnk)$d}else{sigma = svd(xn)$d[1:rnk]})
    irl_noise_scale <- data.table(if(rnk<(p/3)){sigma = irlba(xnprime, nv = rnk)$d}else{sigma = svd(xnprime)$d[1:rnk]})

    return(list(irl_noise, irl_noise_scale))
  }

  irl_noise <- as.data.frame(do.call(cbind,lapply(sq_noise,function(x){x[[1]]})))
  irl_noise_scale <- as.data.frame(do.call(cbind,lapply(sq_noise,function(x){x[[2]]})))

  sigma_n_avg<-rowMeans(irl_noise)
  sigma_n_avg[is.na(sigma_n_avg)]<-0

  sigma_n_avg_scale<-rowMeans(irl_noise_scale)
  sigma_n_avg_scale[is.na(sigma_n_avg_scale)]<-0

  if(scale){
    dat<-data.frame(sigma_n = sigma_n_avg_scale,sigma_a = irl$sigma,dim=1:rnk)
  }else{
    dat<-data.frame(sigma_n = sigma_n_avg,sigma_a = irl$sigma,dim=1:rnk)
  }

  ### Optimization Procedure ###
  gs.model = function(par = par, x = x) {
    b <- par[1]
    c <- par[2]
    d <- par[3]
    sigma_n<-x$sigma_n
    sigma_a<-x$sigma_a
    dim<-x$dim
    gs <- exp(b+c*dim)*(1+d*dim)
    return(gs)
  }

  # run model and compare to true values
  # returns the RMSE
  cost.function <- function(x, par,d) {
    a_obs <- as.vector(x$sigma_a)
    n_pred <- as.vector(gs.model(par = par, x = x))
    dmax <- length(a_obs)
    s <- (a_obs[d:dmax] - n_pred[d:dmax])^2
    n<-length(x$sigma_a[d:dmax])
    w.mean<-NULL
    for (w in 1: n){
      w.mean<-c(w.mean,weighted.mean(s,c(rep(0,w),rep(1,n-w))))
    }
    w.RMSE <- sqrt(mean(w.mean, na.rm = TRUE))
    return(w.RMSE)
  }

  # starting model parameters
  par = c(0,0,0)

  # limits to the parameter space
  # lower <- c(-10,-1000,-100)
  # upper <- c(100,1000,100)

  test = foreach (d =d_start:1,.combine = rbind) %dopar%
  {
    # optimize the model parameters
    optim.par = GenSA(
      par = par,
      d = d,
      fn = cost.function,
      x = dat,
      lower = lower,
      upper = upper,
      control=list(temperature = 4000)
    )
    return(optim.par$par)
  }
  test

  cost = foreach (d =d_start:1,.combine = rbind) %dopar%
  {
    cost.function(par=test[(d_start-d+1),],x=dat,d=d)
  }
  cost
  cost_plot<-ggplot() + theme_minimal() + xlab("Dimension") + ylab("Weighted RMSE") +
    ggtitle("Cost function") +
    geom_line(aes(x = d_start:1, y = cost),colour="black",alpha=0.6) + geom_point(aes(x = d_start:1, y = cost,color = "red"))

  ### Output plots for each fitted model ###
  Plots_To_Save<-list()

  for(d in d_start:1){
    # code for generating each plot is unchanged
    PlottoSave<-ggplot() + theme_minimal() + xlab("Dimension") + ylab("Avg Sigma Noise") +
      theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=12),
            axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),
            plot.title = element_text(hjust = 0.5,face="bold"),legend.position="bottom",
            legend.title = element_text(colour="black", size=12,face="bold"),legend.text=element_text(size=12,face="bold")) +
      ggtitle(paste0("exp(b+c*dim)*(1+d*dim) with dim = ",d,", cost = ", round(cost[(d_start-d+1)],4))) + # ylim(0,250000)+
      geom_line(aes(x = dim, y = sigma_a),dat,colour="black",alpha=0.6) + geom_point(aes(x = dim, y = sigma_a,color = "red"),dat) + geom_text_repel(aes(x = dim, y = sigma_a, label = dim, color = "red"),dat)+
      geom_line(aes(x = dim, y = sigma_n),dat,colour="black",alpha=0.6) + geom_point(aes(x = dim, y = sigma_n,color = "green"),dat) +
      geom_line(aes(x = dat$dim, y = gs.model(par=test[(d_start-d+1),],x=dat)),colour="black",alpha=0.6) + geom_point(aes(x = dat$dim, y = gs.model(par=test[(d_start-d+1),],x=dat),color = "cyan")) +
      scale_color_identity(name = "Model fit",
                           breaks = c("red","green","cyan"),
                           labels = c("sigma_a","sigma_n","exp(b+c*x)-d*x"),
                           guide = "legend") +
      guides(colour = guide_legend(override.aes = list(size=3)))

    Plots_To_Save[[(d_start-d+1)]] <- ggplotGrob(PlottoSave)
  }

  Grid_Plots <- grid.arrange(grobs=Plots_To_Save, ncol=5)
  Grid_Plots


  ### Output Results ###
  result <-list(
      X = X,
      n = n,
      p = p,
      scale = scale,
      d_start = d_start,
      rnk = rnk,
      lower = lower,
      upper = upper
    )
  result$sigma_a <- dat$sigma_a
  result$sigma_n <- dat$sigma_n

  result$parameters <- tibble(b=test[,1],c=test[,2],d=test[,3])
  result$cost <- cost
  result$cost_plot <- cost_plot
  result$Grid_Plots <- Grid_Plots
  result
}
