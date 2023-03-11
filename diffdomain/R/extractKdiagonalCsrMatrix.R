#' make matrix norm
#'
#' @param spsCsrMat  matrix
#'
#' @return matrix
#' @export
#'
#' @examples
#' extractKdiagonalCsrMatrix(matrix)

extractKdiagonalCsrMatrix <- function(spsCsrMat)
{
  #nbin <- nrow(spsCsrMat)
  nonZeroIndex <- NonZero(spsCsrMat)
  nonZeroIndex_row <- nonZeroIndex$row; nonZeroIndex_col <- nonZeroIndex$col

  distkey <- c(); distcolIndex <-c(); distrowIndex <-c(); distvalue <-c()

  for(i in 1:length(nonZeroIndex_row))
  {
    rowIndex <- nonZeroIndex_row[i]
    colIndex <- nonZeroIndex_col[i]
    if(colIndex >= rowIndex)
    {
      distkey <- append(distkey,colIndex - rowIndex)
      distcolIndex <- append(distcolIndex,colIndex)
      distrowIndex <- append(distrowIndex,rowIndex)
      distvalue <- append(distvalue,spsCsrMat[rowIndex,colIndex])
      #print(distvalue)
    }
  }
  UpperMatrixIndex <- data.frame(key = distkey,row = distrowIndex,col = distcolIndex,value = distvalue)
  contactBydistance <- split(UpperMatrixIndex,UpperMatrixIndex$key)
  return(contactBydistance)
}
