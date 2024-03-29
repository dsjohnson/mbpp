#' @title Fit asymptotic approximation of hierarchical N-mixture
#' model for estimation of northern fur seal pup production.
#' @param det_formula formula for the detection model
#' @param avail_formula formula for the avaialbility model
#' @param mark_data Data frame providing pups marked and dead pup counts
#' for each site
#' @param resight_data Data frame providing marked and unmarked resight counts
#' @param par Optional start value specification.
#' @param ... additional arguments passed to \link{R2jags::jags} such as
#' n.chains, n.iter, n.burnin, or n.thin
#' for optimization
#' @author Devin S. Johnson
#' @import dplyr
#' @export

jags_mbpp <- function(dp_formula=~0+log(harem_bulls),
                      avail_formula=~rcode,
                      mark_data, resight_data, par, ...){
  add_args <- list(...)
  data_list = list(
    M=mark_data$M,
    m=resight_data$m,
    n=resight_data$n,
    D=mark_data$deadpups,
    X_a = {
      X <- model.matrix(avail_formula, resight_data)
      idx_avail_occ <- !duplicated(X);
      X <- X[idx_avail_occ,]
      idx <- !(apply(X,2,var)==0 & apply(X,2,mean)!=1)
      X[,idx]
    },
    X_d = {
      X <- model.matrix(dp_formula, mark_data)
      idx <- !(apply(X,2,var)==0 & apply(X,2,mean)!=1)
      X[,idx]
    },
    r_idx = resight_data %>% dplyr::select(rcode, resample) %>% #distinct() %>%
      pull(rcode) %>% {as.integer(factor(.))},
    ro_idx = resight_data %>% {paste0(.$rcode,.$resample)} %>% {as.integer(factor(.))},
    sig_a_idx = attr(model.matrix(avail_formula, resight_data), "assign") %>% as.integer() %>% {.-min(.)+1}
  ) %>% {append(., list(
    n_sig_a = max(.$sig_a_idx)
  ))}

  if(!missing(par)){
    init_list <- par
  } else{
    init_list <- list(
      lambda=as.integer(data_list$M/0.1) - data_list$M,
      tau = rep(0.10, length(data_list$M)),
      theta_a = c(0.5, rep(0, ncol(data_list$X_a)-1)),
      sig_a = c(1, rep(1, max(data_list$sig_a_idx)-1)),
      theta_d = c(0, rep(0, ncol(data_list$X_d)-1))
    )
  }

  par_list <- c(
    "lambda", "tau", "D", "N", "N_tot", "theta_a", "alpha", "Mj", "Uj",
    "theta_d", "sig_a"
  )

  mod_file <- paste0(system.file(package="mbpp"), "/jags/mbpp_jags.bug")
  # mod_file <- paste0(getwd(), "/inst/jags/mbpp_jags_alt.bug")
  if(!("n.chains" %in% names(add_args))){
    init_list <- rep(list(init_list), 3)
  } else{
    init_list <- rep(list(init_list), add_args[["n.chains"]])
  }

  message("Performing MCMC sampling ...")
  fit <- R2jags::jags(data_list, init_list, par_list, mod_file, ...)
  # fit <- R2jags::jags(data_list, init_list, par_list, mod_file)
  message("Creating parameter and abundance summaries ...")
  ## Abundance
  N_data <- mark_data %>% dplyr::select(rcode) %>%
    mutate(
      N = round(apply(fit$BUGSoutput$sims.list$N, 2, median)),
      se_N = round(apply(fit$BUGSoutput$sims.list$N, 2, sd))
    )
  ## availability

  avail_data <- resight_data %>% dplyr::select(rcode, resample) %>%
    distinct() %>%
    mutate(
      alpha=apply(fit$BUGSoutput$sims.list$alpha, 2, median),
      se_alpha = apply(fit$BUGSoutput$sims.list$alpha, 2, sd)
    )
  ## Detection
  det_data <- resight_data %>%
    mutate(
      delta=apply(fit$BUGSoutput$sims.list$delta, 2, median),
      se_delta = apply(fit$BUGSoutput$sims.list$delta, 2, sd)
    ) %>% dplyr::select(rcode, resample, delta, se_delta)
  ## Sigma_alpha
  rnms <- terms(avail_formula) %>%
  {
    rnms <- attr(.,"term.labels")
    if(attr(.,"intercept")) rnms <- c("(Int)", rnms)
    rnms
  }
  sig_a_data <- tibble(
    term = rnms,
    sigma_alpha = apply(fit$BUGSoutput$sims.list$sig_a, 2, median),
    sigma_alpha_se = apply(fit$BUGSoutput$sims.list$sig_a, 2, sd)
  )
  ## Sigma_delta
  rnms <- terms(det_formula) %>%
  {
    rnms <- attr(.,"term.labels")
    if(attr(.,"intercept")) rnms <- c("(Int)", rnms)
    rnms
  }
  sig_d_data <- tibble(
    term = rnms,
    sigma_alpha = apply(fit$BUGSoutput$sims.list$sig_d, 2, median),
    sigma_alpha_se = apply(fit$BUGSoutput$sims.list$sig_d, 2, sd)
  )
  ## Xi
  xi_data <- tibble(
    xi = apply(fit$BUGSoutput$sims.list$xi, 2, median),
    xi_se = apply(fit$BUGSoutput$sims.list$xi, 2, sd)
  )

  return(
    list(
      N = N_data,
      alpha = avail_data,
      sigma_alpha = sig_a_data,
      delta = det_data,
      sigma_delta = sig_d_data,
      xi = xi_data,
      fitting = list(
        summary=fit$BUGSoutput$summary,
        mcmc_samples=fit$BUGSoutput$sims.list,
        last_values = fit$BUGSoutput$last.values
        )
    )
  )

}

get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
