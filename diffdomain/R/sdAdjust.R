#' Return adjust standard deviation
#'
#' @param numberin  number to be calculated
#'
#' @return adjusted standard deviation
#' @export
#'
#' @examples
#' sdAdjust(c(1,2,3,4,5))
sdAdjust <- function(numberin)
{
  n <-length(numberin)
  if(n >=2)
    sdA <- sd(numberin)*sqrt(n-1)*sqrt(1/n)
  else
    sdA <- 0
  return((sdA))
}
