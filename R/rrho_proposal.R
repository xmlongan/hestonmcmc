#' Random generation according to a proposal of \eqn{\rho}
#'
#' Random generation of \eqn{\rho} according to a proposal whose density is
#' defined in [drho_proposal()]. Noting that the generated \eqn{\rho} can not
#' be -1 or 1.
#'
#' Probability function (\eqn{\rho_0^2 \neq 1})
#' \deqn{F(x) =
#'       \begin{cases}
#'         \frac{1}{2}\frac{(x+1)^2}{\rho_0+1}, &-1< x\le \rho_0\\
#'         \frac{1}{2}(\rho_0+1) + \frac{1}{2}\frac{2-x-\rho_0}{1-\rho_0}
#'         (x-\rho_0), &\rho_0 < x < 1
#'       \end{cases}}
#' *Noting that*: the cases \eqn{\rho_0} = 1 or -1 need special treat,
#' * If \eqn{\rho_0=-1},
#' \eqn{F(x) = \frac{1}{2}\frac{2-x+1}{1+1}(x+1),-1<x<1.}
#' * If \eqn{\rho_0=1},
#' \eqn{F(x) = \frac{1}{2}\frac{(x+1)^2}{1+1}, -1<x<1.}
#'
#' Let \eqn{U} is a uniform random over \[0,1\]. We generate X as
#' \eqn{F(X) = U}, i.e., \eqn{X = F^{-1}(U)}
#' \deqn{\begin{matrix}
#'        P(X\le x) &=& P(F^{-1}(U)\le x) \\
#'                  &=& P(U\le F(x)) \\
#'                  &=& F(x)
#'       \end{matrix}}
#' which assures \eqn{X\sim F(x)}. The solution to \eqn{F(X) = U} as following
#' \deqn{X = \begin{cases}
#'       \sqrt{2U(\rho_0+1)} - 1, & U \le \frac{1}{2}(\rho_0+1)\\
#'       1-\sqrt{2(1-U)(1-\rho_0)}, & U > \frac{1}{2}(\rho_0+1)
#' \end{cases}}
#' Above generator for \eqn{X} works for all cases, i.e., \eqn{\rho_0=-1}, or
#' \eqn{\rho_0=1}, or \eqn{\rho_0\neq -1}, or \eqn{\rho_0\neq 1}. Note that
#' \eqn{X} can be neither 1 nor -1 which is assured by the fact `runif(1,0,1)`
#'  will not generate either 0 or 1.
#'
#' @param rho0 center point of the proposal distribution.
#'
#' @return a scale value.
#' @export
#'
#' @examples
#' rrho_proposal(0.1)
#' rrho_proposal(0.5)
rrho_proposal <- function(rho0) {
  U = stats::runif(1)
  if (U < (rho0+1)/2) {
    rho = sqrt(2*U*(rho0+1)) - 1
  } else {
    rho = 1 - sqrt(2*(1-U)*(1-rho0))
  }
  return(rho)
}
