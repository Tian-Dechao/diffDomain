#' Return matrix from hic file
#'
#' @param chrn  chosen chromosome
#' @param start start bin of chromosome
#' @param end end bin of chromosome
#' @param reso chosen resoltion
#' @param fhic chosen hic file
#' @param hicnorm default is "KR"
#'
#' @return matrix
#' @export
#'
#' @examples
#' contact_matrix_from_hic(1,1060000,1080000,10000,'/lanec1_home/data/GSE63525_cell2014/K562_combined.hic','KR')
contact_matrix_from_hic <- function(chrn,start,end,reso,fhic,hicnorm)
{
  options(scipen = 200)
  options(digits = 8)
  # handel some regions cannot be loaded by straw
  # find the bins for a domain
  domwin <- makewindow2(start,end,reso)
  # create index for domwin
  nb <- length(domwin)
  nbs <- 1:nb
  domwin_dict <- data.frame(domwin = domwin,nb = nbs)
  # create empty matrix
  mat <- matrix(data = NA,nrow=nb,ncol=nb)
  # find the edgelist from .hic
  region <- stringr::str_c(chrn,domwin[1],domwin[nb],sep = ":")
  print(region)

  # straw
  if(substring(fhic,stringr::str_length(fhic)-3,stringr::str_length(fhic))  == '.hic')
  {
    el <- strawr::straw(hicnorm,fhic,region,region,'BP',reso)
    for(i in 1:nrow(el))
    {
      bin0 <- el[,1][i]
      bin1 <- el[,2][i]
      k <- domwin_dict[which(domwin_dict$domwin == bin0),]$nb
      l <- domwin_dict[which(domwin_dict$domwin == bin1),]$nb
      if(k == l)
      {
        mat[k,l] <- el[,3][i]
      }
      else
      {
        mat[k,l] <- el[,3][i]
        mat[l,k] <- el[,3][i]
      }
    }
    return(mat)
  }

  else
  {
    sprintf("Sorry, %s dosn't exist" ,fhic)
    mat = NULL
    return(mat)
  }

}

