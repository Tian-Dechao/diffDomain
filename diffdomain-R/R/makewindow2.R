
#' Return windows between start and end,based on the resolution
#'
#' @param start start bin of chosen chromosome
#' @param end end bin of chosen chromosome
#' @param reso chosen resolution
#'
#' @return wins
#' @export
#'
#' @examples
#' makewindow2(1060000,1080000,10000)
makewindow2 <- function(start,end,reso)
{
  wins = c()

  k1 <- start %/% reso
  k2 <- end %/% reso
  l2 <- end %% reso
  if(l2 >= 0)
    k2 = k2 + 1

  for(i in k1:(k2-1))
  {
    wins <- append(wins,i*reso)
  }
  return(wins)
}
