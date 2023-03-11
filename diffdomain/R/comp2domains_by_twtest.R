#' compare single pair of domains
#'
#' @param chrn chromosomee to compare
#' @param start  start bin of chrn
#' @param end  end bin of chrn
#' @param fhic0  the first hic file to be compared
#' @param fhic1  the  second hic file to be compared
#' @param ofile  the file path of the result
#' @param min_nbin default is 10
#' @param hicnorm  default is "KR"
#' @param reso  resolution to choose,default is 100000
#' @param f  default is 0.5
#'
#' @return result of comparition
#' @export
#'
#' @examples
#' comp2domains_by_twtest('/lanec1_home/data/GSE63525_cell2014/GSE63525_GM12878_insitu_primary_replicate_combined.hic',
#' '/lanec1_home/data/GSE63525_cell2014/GSE63525_K562_combined.hic',
#' ofile = '/lanec2_home/zhangx/DiffCompare/pythonR_compare_result/stdout_R')

comp2domains_by_twtest <- function(chrn,start,end,fhic0,fhic1,
                                   ofile = paste(getwd(),'/stdout.txt',sep = ''),min_nbin = 10,hicnorm = 'KR',reso = 100000,f = 0.5)
{
  options(scipen = 200)
  options(digits = 8)
  mat0 <- contact_matrix_from_hic(chrn,start,end,reso,fhic0,hicnorm)
  mat1 <- contact_matrix_from_hic(chrn,start,end,reso,fhic1,hicnorm)

  if(! is.null(mat0) & ! is.null(mat1))
  {
    # move rows that have more than half NA
    nbins <- compute_nbins(start,end,reso)
    ind0 <-c() ; ind1 <- c()
    for(i in 1:nrow(mat0))
    {
      ind0[i] <- sum(is.na(mat0[,i]))
      ind1[i] <- sum(is.na(mat1[,i]))
    }
    ind0 <- ind0 < nbins * f
    ind1 <- ind1 < nbins * f
    ind = ind0 & ind1

    if(sum(ind) >= min_nbin)
    {
      mat0rmna <- mat0[ind,]
      mat0rmna <- mat0rmna[,ind]
      mat1rmna <- mat1[ind,]
      mat1rmna <- mat1rmna[,ind]

      # compute th difference matrix
      Diffmat <- mat0rmna / mat1rmna
      ind1 = c()
      # print Diffmat
      for(i in 1:nrow(Diffmat))
        ind1[i] <- sum(is.na(Diffmat[,i]))

      Diffmatnorm <- normDiffbyMeanSD(D = Diffmat)

      result <- twtest_formula(Diffmatnorm)

      domname <- paste(sprintf('chr%s',chrn),sprintf(':%s',start),sprintf('-%s',end),sep = '')
      result = data.frame(chrn = chrn,start = start,end = end, region = domname,
                          lambdan = result$lambdan, p = result$p,nd = result$nd)
      #write.table(data.frame(result),file = paste(path,'/stdout.txt',sep = ''),sep = '\t')
    }
    else
    {
      domname <- paste(sprintf('chr%s',chrn),sprintf(':%s',start),sprintf('-%s',end),sep = '')
      result = data.frame(chrn = chrn,start = start,end = end,region = domname,lambdan = NA,p =  NA, nd = NA)
      print('The length of this TAD is too small at this resolution to be calculated !')
    }
  }
  else
  {
    domname <- paste(sprintf('chr%s',chrn),sprintf(':%s',start),sprintf('-%s',end),sep = '')
    print('The matrix is too sparse at this resolution to be calculated !')

    result = data.frame(chrn = chrn,start = start,end = end,region = domname,
                        lambdan = NA, p =  NA, nd = NA)
  }

  write.table(data.frame(result),ofile,sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)
  return(result)
}
