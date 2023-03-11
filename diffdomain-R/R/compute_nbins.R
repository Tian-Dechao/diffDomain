#' return how many bins between start and end, based on the chosen resolution
#'
#' @param start start bin of chromosome
#' @param end end bin of chromosome
#' @param reso chosen resolution
#'
#' @return number of bins
#' @export
#'
#' @examples
#' compute_nbins(1060000,1080000,10000)

compute_nbins <- function(start,end,reso)
{
  domwin <- makewindow2(start,end,reso)
  # create index for domwin
  nb <- length(domwin)
  return(nb)
}
