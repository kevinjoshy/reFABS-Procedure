library(metap)
library(BLMA)

reFABSp <- function(pvals)
{
  res <- numeric(length = length(pvals))
  for(i in 1:length(pvals))
  {
    if(length(pvals[i][[1]]) == 1)
    {
      res[i] <- pvals[i][[1]][1]
    }
    if(length(pvals[i][[1]]) > 1)
    {
      z = allmetap(pvals[[i]], method = c('logitp', 'maximump', 'minimump', 'sumlog', 'sump', 'sumz'))
      minz = min(as.numeric(z$p))
      x = fisherMethod(pvals[[i]])
      res[i] = min(minz, x)
    }
  }
  return(res)
}



