#' Find those whose value equal zero
#'
#' @param D matrix
#'
#' @return whose value equals zero and their rows and cols
#' @export
#'
#' @examples
#'NonZero(input)
NonZero <- function(D)
{
  nonZeroIndex_row <- c(); nonZeroIndex_col <- c()
  for(i in 1:nrow(D))
  {
    for(j in 1:nrow(D))
    {
      if(D[i,j] != 0 | is.na(D[i,j]) | D[i,j] == Inf)
      {
        nonZeroIndex_row <- append(nonZeroIndex_row,i)
        nonZeroIndex_col <- append(nonZeroIndex_col,j)
      }
    }
    result <- list(row = nonZeroIndex_row,col = nonZeroIndex_col)
  }
  return(result)
}
