#' Return twtest result of matrix
#'
#' @param Mat  input matrix
#'
#' @return test result of chosen matrix
#' @export
#'
#' @examples
#' twtest_formula(matrix
twtest_formula <- function(Mat)
{
  options(scipen = 200)
  options(digits = 8)

  nd <- nrow(Mat)
  Mat <- Mat /sqrt(nd)
  tryCatch(
    {
      w <- eigen(Mat)$values ; v <- eigen(Mat)$vectors
      lambdan <- max(abs(w))
      lambdan <- (lambdan-2) * nd^(2/3)
      twpercent <- RMTstat::ptw(q = round(lambdan,digits = 5),beta = 1)
      p = 1- twpercent
      result <- list(nd = nd,lambdan = round(lambdan,digits = 5),p = round(p,digits = 7))
      return(result)
    },
    error=function()
    {
      message('An Error Occurred')
      print(e)
      np = NA; lambdan = NA; p = NA
      result <- list(nd = nd,lambdan = lambdan,p = p)
      return(result)
    },
    warning=function()
    {
      message('A Warning Occurred')
      print(w)
      np = NA; lambdan = NA; p = NA
      result <- list(nd = nd,lambdan = lambdan,p = p)
      return(result)
    }
  )

}
