#' compare multiple domains
#'
#' @param fhic0 the first hic file to be compared
#' @param fhic1 the  second hic file to be compared
#' @param tadlist_of_hic0.bed domian list
#' @param ofile the file path of the result
#' @param min_nbin default is 10
#' @param hicnorm default is "KR"
#' @param reso  resolution to choose,default is 100000
#' @param f default is 0.5
#'
#' @return result of single pair comparition
#' @export
#'
#' @examples
#'comp2domains_by_twtest_parallel(
#''/lanec1_home/data/GSE63525_cell2014/GSE63525_GM12878_insitu_primary_replicate_combined.hic',
#''/lanec1_home/data/GSE63525_cell2014/GSE63525_KBM7_combined.hic',
#''/lanec1_home/data/tads/ArrowheadDomainlist/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt',
#'ofile = '/lanec2_home/zhangx/DiffCompare/pythonR_compare_result/GM12878_KBM7_stdout_R')

comp2domains_by_twtest_parallel <- function(fhic0,fhic1,tadlist_of_hic0.bed,
                                            ofile=paste(getwd(),'/stdout',sep = ''),min_nbin = 10,hicnorm = 'KR',reso = 100000,f = 0.5)
{
  options(scipen = 200)
  options(digits = 8)

  tadlist <- read.table(tadlist_of_hic0.bed,header = T)
  result = data.frame()
  for(i in 1:nrow(tadlist))
  {
    result = rbind(result,comp2domains_by_twtest(chrn = tadlist[i,][,1],
                                                 start = tadlist[i,][,2],
                                                 end = tadlist[i,][,3],
                                                 fhic0,fhic1,
                                                 min_nbin = 10,hicnorm = 'KR',reso = 100000,f = 0.5))


  }
  result <- na.omit(result)
  write.table(result, ofile,sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)
  return(result)
}
