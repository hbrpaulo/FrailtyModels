#' Flexible frailty model utilities
#'
#' Helper functions for fitting a frailty model with customizable baseline
#' hazards and empirical Bayes estimation.
#'
#' @examples
#' \dontrun{
#' # Script depends on external data; see functions below for usage
#' }
NULL

# Flexible frailty model with customizable baseline hazard (H0) and EB frailty estimation (with covariate integration)

dist_chosen <- "weibull"  # Options: "weibull", "exponential", "gompertz", "gamma", "lognormal", "loglogistic", "piecewise_exp"
breaks_pw <- c(0, 50, 100, 200)

data_real <- read.csv("~/GitHub/pipges/petrobras/ICV apresentacao 02/dados/clusters_2_3.csv")

surv_time <- data_real$tt_fi_dias
event <- as.numeric(!data_real$censura)
covariates <- data_real %>% select(tt_fi_dias, fabr)
X_real <- model.matrix(tt_fi_dias~., covariates)

baseline_functions <- list(
  weibull = list(
    H0 = function(t, shape, scale) (t/scale)^shape,
    h0 = function(t, shape, scale) (shape/scale^shape) * t^(shape-1)
  ),
  exponential = list(
    H0 = function(t, rate) rate * t,
    h0 = function(t, rate) rep(rate, length(t))
  ),
  gompertz = list(
    H0 = function(t, a, b) (a/b) * (exp(b*t)-1),
    h0 = function(t, a, b) a * exp(b*t)
  ),
  gamma = list(
    H0 = function(t, shape, rate) -log(1 - pgamma(t, shape=shape, rate=rate)),
    h0 = function(t, shape, rate) dgamma(t, shape=shape, rate=rate)/(1-pgamma(t,shape=shape,rate=rate))
  ),
  lognormal = list(
    H0 = function(t, meanlog, sdlog) -log(1 - plnorm(t, meanlog=meanlog, sdlog=sdlog)),
    h0 = function(t, meanlog, sdlog) dlnorm(t, meanlog=meanlog, sdlog=sdlog)/(1-plnorm(t,meanlog=meanlog,sdlog=sdlog))
  ),
  loglogistic = list(
    H0 = function(t, shape, scale) -log(1/(1+(t/scale)^shape)),
    h0 = function(t, shape, scale) ((shape/scale)*(t/scale)^(shape-1))/(1+(t/scale)^shape)
  ),
  piecewise_exp = list(
    H0 = function(t, breaks, rates) {
      sapply(t, function(x){
        idx <- findInterval(x, breaks, rightmost.closed=TRUE)
        cumhaz <- sum(diff(c(0,breaks[1:idx]))*rates[1:idx])
        if(idx < length(rates)) cumhaz <- cumhaz+(x-breaks[idx])*rates[idx+1]
        cumhaz
      })
    },
    h0 = function(t, breaks, rates) {
      sapply(t, function(x){ rates[findInterval(x,breaks,rightmost.closed=TRUE)+1] })
    }
  )
)

#' Log-likelihood for flexible frailty model
#'
#' Computes the negative log-likelihood for a frailty model with
#' distribution-specific baselines and covariate effects.
#'
#' @param params Parameter vector containing baseline, frailty variance,
#'   and covariate effects.
#' @param times Follow-up times.
#' @param event Event indicators.
#' @param X Covariate matrix.
#' @param dist Baseline distribution name.
#' @param breaks Break points for piecewise exponential baseline.
#'
#' @return Numeric scalar negative log-likelihood.
#' @examples
#' \dontrun{
#' frailty_flexible_loglik(rep(0.1, 5), surv_time, event, X_real)
#' }
frailty_flexible_loglik <- function(params, times, event, X, dist=dist_chosen, breaks=breaks_pw){
  frailty_var <- params[3]
  p <- ncol(X)
  beta_shape <- params[4:(3+p)]
  beta_scale <- params[(4+p):(3+2*p)]
  alpha_gamma <- 1/frailty_var
  
  if(dist=="weibull"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    scale_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$weibull$H0(times, shape_ind, scale_ind)
    h0 <- baseline_functions$weibull$h0(times, shape_ind, scale_ind)
  } else if(dist=="exponential"){
    rate_ind <- exp(log(params[1])+X%*%beta_scale)
    H0 <- baseline_functions$exponential$H0(times, rate_ind)
    h0 <- baseline_functions$exponential$h0(times, rate_ind)
  } else if(dist=="gompertz"){
    a_ind <- exp(log(params[1])+X%*%beta_shape)
    b_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$gompertz$H0(times, a_ind, b_ind)
    h0 <- baseline_functions$gompertz$h0(times, a_ind, b_ind)
  } else if(dist=="gamma"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    rate_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$gamma$H0(times, shape_ind, rate_ind)
    h0 <- baseline_functions$gamma$h0(times, shape_ind, rate_ind)
  } else if(dist=="lognormal"){
    meanlog_ind <- log(params[1])+X%*%beta_shape
    sdlog_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$lognormal$H0(times, meanlog_ind, sdlog_ind)
    h0 <- baseline_functions$lognormal$h0(times, meanlog_ind, sdlog_ind)
  } else if(dist=="loglogistic"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    scale_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$loglogistic$H0(times, shape_ind, scale_ind)
    h0 <- baseline_functions$loglogistic$h0(times, shape_ind, scale_ind)
  } else if(dist=="piecewise_exp"){
    rates <- exp(log(params[1:(length(breaks)-1)])+X%*%beta_scale[1:(length(breaks)-1)])
    H0 <- baseline_functions$piecewise_exp$H0(times, breaks, rates)
    h0 <- baseline_functions$piecewise_exp$h0(times, breaks, rates)
  } else stop("Unsupported baseline")
  
  S <- (1+H0/alpha_gamma)^(-alpha_gamma)
  h <- h0/(1+H0/alpha_gamma)
  if(any(S<=0)||any(event==1 & h<=0)) return(Inf)
  -sum(event*(log(h)+log(S))+(1-event)*log(S))
}

#' Empirical Bayes frailty estimate
#'
#' Computes the posterior mean frailty for each observation under the
#' specified baseline hazard and covariate effects.
#'
#' @inheritParams frailty_flexible_loglik
#'
#' @return Numeric vector of expected frailty values.
#' @examples
#' \dontrun{
#' EB_frailty_flexible(rep(0.1, 5), surv_time, event, X_real)
#' }
EB_frailty_flexible <- function(params, times, event, X, dist=dist_chosen, breaks=breaks_pw){
  frailty_var <- params[3]
  alpha_gamma <- 1/frailty_var
  p <- ncol(X)
  beta_shape <- params[4:(3+p)]
  beta_scale <- params[(4+p):(3+2*p)]
  
  if(dist=="weibull"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    scale_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$weibull$H0(times, shape_ind, scale_ind)
  } else if(dist=="exponential"){
    rate_ind <- exp(log(params[1])+X%*%beta_scale)
    H0 <- baseline_functions$exponential$H0(times, rate_ind)
  } else if(dist=="gompertz"){
    a_ind <- exp(log(params[1])+X%*%beta_shape)
    b_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$gompertz$H0(times, a_ind, b_ind)
  } else if(dist=="gamma"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    rate_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$gamma$H0(times, shape_ind, rate_ind)
  } else if(dist=="lognormal"){
    meanlog_ind <- log(params[1])+X%*%beta_shape
    sdlog_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$lognormal$H0(times, meanlog_ind, sdlog_ind)
  } else if(dist=="loglogistic"){
    shape_ind <- exp(log(params[1])+X%*%beta_shape)
    scale_ind <- exp(log(params[2])+X%*%beta_scale)
    H0 <- baseline_functions$loglogistic$H0(times, shape_ind, scale_ind)
  } else if(dist=="piecewise_exp"){
    rates <- exp(log(params[1:(length(breaks)-1)])+X%*%beta_scale[1:(length(breaks)-1)])
    H0 <- baseline_functions$piecewise_exp$H0(times, breaks, rates)
  } else stop("Unsupported baseline")
  
  (alpha_gamma+event)/(alpha_gamma+H0)
}

if(dist_chosen=="piecewise_exp"){
  init_params <- c(rep(0.01,length(breaks_pw)-1),0.3,rep(0,ncol(X_real)),rep(0,ncol(X_real)))
} else {
  init_params <- c(2,1.5,0.3,rep(0,ncol(X_real)),rep(0,ncol(X_real)))
}

opt_real <- optim(init_params,
                  frailty_flexible_loglik,
                  times=surv_time,event=event,
                  X=X_real,dist=dist_chosen,
                  breaks=breaks_pw,
                  method="L-BFGS-B",
                  lower=rep(0.001,length(init_params)),
                  hessian=TRUE)

est_real <- opt_real$par
vcov_real <- solve(opt_real$hessian)
se_real <- sqrt(diag(vcov_real))
z_val <- qnorm(1-0.05/2)
lower_real <- est_real - z_val*se_real
upper_real <- est_real + z_val*se_real

table_real <- data.frame(
  Parameter=c("Baseline Params","Frailty Var",paste0("Beta_Shape",1:ncol(X_real)),
              paste0("Beta_Scale", 1:ncol(X_real))), Estimate=est_real,
  Lower=lower_real, Upper=upper_real, SE=se_real)
print(table_real)

EB_frailty_real <- EB_frailty_flexible(est_real,surv_time,event,X_real,dist=dist_chosen,breaks=breaks_pw)
ggplot(data.frame(frailty=EB_frailty_real),aes(x=frailty))+
  geom_histogram(aes(y=after_stat(density)),bins=25,fill="skyblue",color="darkblue",alpha=0.6)+
  geom_density(linetype="dashed",color="darkgreen",linewidth=1.2)+
  labs(title="Empirical Bayes Frailty Estimates - Real Data",x="Frailty (Z)",y="Density",subtitle=paste("Distribution:",str_to_title(dist_chosen)))+
  theme_minimal()
