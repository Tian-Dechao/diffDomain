#' Return nomalization matrix based on its mean and standard deviation
#'
#' @param D  matrix
#'
#' @return normalized matrix
#' @export
#'
#' @examples
#' normDiffbyMeanSD(matrix)
normDiffbyMeanSD <- function(D)
{
  options(scipen = 200)
  options(digits = 8)

  D <- log(D)
  contactByDistanceDiff <- extractKdiagonalCsrMatrix(D) # output delete inf

  # imputation of nan and inf
  # part0: get the median and maximum for each off-diagonal
  key <- c();
  contactMedian <- data.frame(); contactMax <- data.frame(); contactMean <- data.frame(); contactStd <- data.frame(

  )
  for(k in 1:length(contactByDistanceDiff))
  {
    val <- contactByDistanceDiff[[k]]$value
    indnan <- is.na(val)
    indinf <- is.infinite((val))

    if(sum(indnan) > 0 & sum(indinf) > 0)
    {
      ind <- indnan | indinf
      val1 <- val[!ind] # delete NA and Inf
      DataFramesum <- data.frame(key = k,contactmedian = median(val1))
      contactMedian <- rbind(contactMedian,DataFramesum)

      DataFramesum <- data.frame(key = k,contactmean = mean(val1))
      contactMean <- rbind(contactMean,DataFramesum)

      DataFramesum <- data.frame(key = k,contactmax = max(val1))
      contactMax <- rbind(contactMax,DataFramesum)

      DataFramesum <- data.frame(key = k,contactstd = sdAdjust(val1))
      contactStd <- rbind(contactStd,DataFramesum)
    }
    else if(sum(indnan) > 0 & sum(indinf) == 0)
    {
      val1 = val[!indnan] # delete NA

      DataFramesum <- data.frame(key = k,contactmedian = median(val1))
      contactMedian <- rbind(contactMedian,DataFramesum)

      DataFramesum <- data.frame(key = k,contactmean = mean(val1))
      contactMean <- rbind(contactMean,DataFramesum)

      DataFramesum <- data.frame(key = k,contactstd = sdAdjust(val1))
      contactStd <- rbind(contactStd,DataFramesum)
    }

    else if (sum(indnan) == 0 & sum(indinf) > 0)
    {
      val1 = val[!indinf] # delete Inf
      DataFramesum <- data.frame(key = k,contactmean = mean(val1))
      contactMean <- rbind(contactMean,DataFramesum)

      DataFramesum <- data.frame(key = k,contactmax = max(val1))
      contactMax <- rbind(contactMax,DataFramesum)

      DataFramesum <- data.frame(key = k,contactstd = sdAdjust(val1))
      contactStd <- rbind(contactStd,DataFramesum)
    }
    else
    {
      DataFramesum <- data.frame(key = k,contactmean = mean(val))
      contactMean <- rbind(contactMean,DataFramesum)

      DataFramesum <- data.frame(key = k,contactstd = sdAdjust(val))
      contactStd <- rbind(contactStd,DataFramesum)
    }
  }

  indnan <- is.na(D)
  indnan[indnan] <- 1
  indnan <- NonZero(indnan)
  indr <- indnan$row; indc <- indnan$col
  for(i in 1:length(indr))
  {
    k <- abs(indr[i]-indc[i]) + 1
    D[indr[i],indc[i]] <- contactMedian[which(contactMedian$key == k),]$contactmedian
  }

  indinf <- is.infinite(D)
  indinf[indinf] <- 1
  indinf <- NonZero(indinf)
  indr <- indinf$row; indc <- indinf$col
  for(i in 1:length(indr))
  {
    k <- abs(indr[i]-indc[i]) + 1
    D[indr[i],indc[i]] <- contactMax[which(contactMax$key == k),]$contactmax
  }

  # sbstract mean and dividing by sd
  for(i in 1:nrow(contactStd))
  {
    if(is.na(contactStd$contactstd[i]) | contactStd$contactstd[i] == 0)
      contactStd$contactstd[i] <- 1
  }

  for(i in 1:nrow(D))
  {
    for(j in 1:nrow(D))
    {
      k <- abs(i-j) + 1
      D[i,j] <- (D[i,j]-contactMean$contactmean[k]) / contactStd$contactstd[k]
    }
  }

  return(D)
}
