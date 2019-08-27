#' @title Fit asymptotic approximation of hierarchical N-mixture
#' model for estimation of northern fur seal pup production.
#' @param det_formula formula for the detection model
#' @param avail_formula formula for the avaialbility model
#' @param mark_data Data frame providing pups marked and dead pup counts
#' for each site
#' @param resight_data Data frame providing marked and unmarked resight counts
#' @param par Optional start value specification.
#' @param ... additional arguments passed to nlminb() or TMB::MakeADFun()
#' for optimization
#' @author Devin S. Johnson
#' @import TMB
#' @importFrom numDeriv grad
#' @export

asym_mbpp <- function(
  det_formula=~rcode*resample*observer,
  avail_formula=~rcode,
  mark_data, resight_data, par, ...){

  data_list = list(
    omega=c(
      factor(resight_data$rcode) %>% as.integer(),
      factor(resight_data$resample) %>% as.integer(),
      factor(resight_data$observer) %>% as.integer()
    ) %>% matrix(ncol=3)-1,
    M=mark_data$M,
    m=resight_data$m,
    u=resight_data$u,
    resight_data$u,
    D=mark_data$deadpups,
    X_a = {
      X <- model.matrix(avail_formula, resight_data)
      idx <- !(apply(X,2,var)==0 & apply(X,2,mean)!=1)
      X[,idx]
    },
    X_d = {
      X <- model.matrix(det_formula, resight_data)
      idx <- !(apply(X,2,var)==0 & apply(X,2,mean)!=1)
      X[,idx]
    },
    sig_a_idx = attr(model.matrix(avail_formula, resight_data), "assign") %>% as.integer() %>% {.-min(.)},
    sig_d_idx = attr(model.matrix(det_formula, resight_data), "assign") %>% as.integer() %>% {.-min(.)} #int
  )

  if(missing(par)){
    par_list <- list(
      log_lambda=log(data_list$M/0.1),
      logit_xi = log(0.02),
      logit_tau = rep(qlogis(0.10), length(data_list$M)),
      theta_a = rep(0, ncol(data_list$X_a)),
      log_sig_a = c(log(5), rep(log(1), max(data_list$sig_a_idx))),
      theta_d = rep(0, ncol(data_list$X_d)),
      log_sig_d = c(log(5),rep(log(1), max(data_list$sig_d_idx)))
    )
  } else{
    par_list <- par
  }
  ## Load tmb ddl
  pth <- paste0(system.file(package="mbpp"), "/tmb/")
  chk <- as.logical(suppressMessages(mbpp:::check_compile()))
  if(chk){
    message("TMB source code must be compiled...")
    TMB::compile(paste0(pth, "mbpp.cpp"))
  }
  dyn.load(TMB::dynlib(paste0(pth, "mbpp")))

  obj <- TMB::MakeADFun(data_list, par_list,
                        random=c("theta_a","theta_d"),
                        DLL="mbpp", ...)
  message("Optimizing integrated likelihood ...")
  opt <- nlminb(obj$par, obj$fn, obj$gr, ...)
  if(opt$convergence!=0){
    message("Optimization did not converge!")
    return(list(obj=obj, opt=opt, data_list=data_list, par_list=par_list))
  } else{
    message("Creating parameter and abundance summaries ...")
    sdrep <- TMB::sdreport(obj,getJointPrecision=TRUE)
    Cmat <- sdrep$jointPrecision %>% solve()
    summ_sdrep <- summary(sdrep)
    nms <- rownames(summ_sdrep)
    nms_c <- rownames(Cmat)
    summ_sdrep <- as.data.frame(summ_sdrep)

    ## Abundance
    V_miss <- summ_sdrep %>% .[(nms=="miss"),"Std. Error"]
    E_miss <- summ_sdrep %>% .[(nms=="miss"),"Estimate"]
    N_data <- mark_data %>% dplyr::select(rcode) %>%
      mutate(
        N = round(summ_sdrep %>% .[(nms=="N_est"),"Estimate"]),
        se_N = round(sqrt(V_miss^2 + E_miss))
      )
    ## availability
    X_a <- data_list$X_a
    v_alpha <- {X_a%*% Cmat[nms_c=="theta_a",nms_c=="theta_a"] %*%t(X_a)} %>%
      diag()
    logit_alpha <- {X_a%*%summ_sdrep[nms=="theta_a",1]}
    gr_alpha <- numDeriv::grad(plogis, logit_alpha)
    v_alpha <- gr_alpha^2 * v_alpha
    avail_data <- resight_data %>%
      mutate(
        alpha=plogis(logit_alpha) %>% as.vector(),
        se_alpha = sqrt(v_alpha)
      ) %>% dplyr::select(rcode, resample, alpha, se_alpha) %>%
      distinct()
    ## Detection
    X_d <- data_list$X_d
    v_delta <- {X_d%*% Cmat[nms_c=="theta_d",nms_c=="theta_d"] %*%t(X_d)} %>%
      diag()
    logit_delta <- {X_d%*%summ_sdrep[nms=="theta_d",1]}
    gr_delta <- numDeriv::grad(plogis, logit_delta)
    v_delta <- gr_delta^2 * v_delta
    det_data <- resight_data %>%
      mutate(
        delta=plogis(logit_delta) %>% as.vector(),
        se_delta = sqrt(v_delta)
      ) %>% dplyr::select(rcode, resample, delta, se_delta)
    ## Sigma_alpha
    rnms <- terms(avail_formula) %>%
    {
      rnms <- attr(.,"term.labels")
      if(attr(.,"intercept")) rnms <- c("(Int)", rnms)
      rnms
    }
    sig_a_data <- summ_sdrep[nms=="sig_a",] %>%
      `colnames<-`(c("sigma_alpha","sigma_alpha_se")) %>%
      `rownames<-`(NULL) %>% as_tibble() %>%
      mutate(term = rnms) %>% .[,c(3,1,2)]

    ## Sigma_delta
    rnms <- terms(det_formula) %>%
      {
        rnms <- attr(.,"term.labels")
        if(attr(.,"intercept")) rnms <- c("(Int)", rnms)
        rnms
      }
    sig_d_data <- summ_sdrep[nms=="sig_d",] %>%
      `colnames<-`(c("sigma_delta","sigma_delta_se")) %>%
      `rownames<-`(NULL) %>% as_tibble() %>%
      mutate(term = rnms) %>% .[,c(3,1,2)]
    ## xi
    xi_data <- summ_sdrep %>% filter(nms=="xi") %>%
      rename(xi=Estimate, xi_se=`Std. Error`)
    ## Raw par
    order <- nms_c %>% unique()
    par = NULL
    for(i in 1:length(order)){
      tmp <- summ_sdrep %>% filter(nms==order[i]) %>% pull(Estimate)
      par <- c(par, tmp)
    }
    names(par) <- nms_c
    return(
      list(
        N = N_data,
        alpha = avail_data,
        sigma_alpha = sig_a_data,
        delta = det_data,
        sigma_delta = sig_d_data,
        xi = xi_data,
        data_list = data_list,
        fitting = list(obj=obj, opt=opt, par=par, covmat=Cmat)
      )
    )
  }



}
