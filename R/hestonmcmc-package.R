#' @keywords internal
#' @author Yanfeng Wu <yanfengwu13@fudan.edu.cn>
#' 
#' @details 
#' We try to implement the MCMC as close as to that in Eraker et al. (2003) and
#' Johannes and Polson (2010). However, there are some details are implemented
#' according to my judgment, because of the vague (maybe errors) in the two 
#' papers.
#' 
#' @references
#' Eraker, B., Johannes, M., & Polson, N. (2003). The impact of jumps in 
#' volatility and returns. *The Journal of Finance*, 58(3), 1269-1300.
#' 
#' Johannes, M., & Polson, N. (2010). MCMC methods for continuous-time 
#' financial econometrics. In 
#' *Handbook of Financial Econometrics: Applications* (pp. 1-72). Elsevier.
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib hestonmcmc, .registration = TRUE
## usethis namespace: end
NULL
