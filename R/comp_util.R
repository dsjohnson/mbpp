#' @title Compile package TMB function
#' @param model Details which model to compile
#' @details This function is only used once after the
#' package is installed to compile the c++ TMB code for
#' asymptotic abundance estimation. Simply call the function
#' to compile the source code.
#' @importFrom TMB compile
#' @export
compile_mbpp_tmb <- function(){
  ftmb <- system.file("tmb", "mbpp.cpp", package="mbpp")
  TMB::compile(ftmb)
}

check_compile <- function(){
  path <- paste0(system.file(package="mbpp"), "/tmb/")
  fso <- paste0(path, "mbpp.so")
  fdll <- paste0(path, "mbpp.dll")
  if(file.exists(fso) | file.exists(fdll)){
    return(0)
  } else{
    message("\nLooks like the TMB code has not been compiled.\nPlease run 'compile_mbpp_tmb()'\n")
    return(1)
  }
}
