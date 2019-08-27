#' @title Parameteric bootstrap sample for abundance
#' @details This function uses the asymtotic model fit from TMB to generate
#' a bootstrap sample for parameter inference
#' @param object A fitted model object from aysm_mbpp
#' @param size Number of bootstrap samples to draw
#' @importFrom dplyr select
#' @export
#'
boot_asym <- function(object, size=10000){
  par <- object$fitting$par
  Cmat <- object$fitting$covmat
  dl <- object$data_list
  idx <- colnames(Cmat) %in% c("log_lambda","logit_xi","logit_tau")
  Cmat <- Cmat[idx,idx]
  par <- par[idx]
  boot <- mvtnorm::rmvnorm(size, par, Cmat)
  lam <- boot[,colnames(boot)%in%c("log_lambda")] %>% exp()
  tau <- boot[,colnames(boot)%in%c("logit_tau")] %>% plogis()
  xi <- boot[,colnames(boot)%in%c("logit_xi")] %>% plogis()
  lam_U <- (1-xi)*(1-tau)*lam
  U <- rpois(ncol(lam_U)*size, lam_U) %>% matrix(nrow=size)
  lam_D = xi*lam
  D <- rpois(ncol(lam_D)*size, lam_D) %>% matrix(nrow=size)
  Nmat <- U
  for(i in 1:ncol(U)){
    Nmat[,i] <- Nmat[,i] + ifelse(is.na(dl$D[i]), 1, 0)*D[,i] +
      ifelse(!is.na(dl$D[i]), dl$D[i], 0)+ dl$M[i]
  }
  colnames(Nmat) <- object$N$rcode
  N_data <- object$N %>% dplyr::select(rcode) %>%
    mutate(
      N = round(apply(Nmat, 2, median)),
      se_N = round(apply(Nmat, 2, sd))
    )
  return(
    list(
      N = N_data,
      N_boot = Nmat
    )
  )
}
