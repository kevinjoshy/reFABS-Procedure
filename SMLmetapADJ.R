library(metap)
# ignore meanz, meanp, na, NaN
reFABSp <- function(pvals) {
  z = allmetap(pvals, method = c('logitp', 'maximump', 'minimump', 'sumlog', 'sump', 'sumz'))
  min(as.numeric(z$p))
}



