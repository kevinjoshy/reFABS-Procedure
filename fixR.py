import rpy2.robjects as robjects
import os

os.environ['R_HOME'] = '/Users/joshyk/Library/R/4.0/library'
# Combines the three different p-vals together
def reFABScalc(pvals):
    r = robjects.r
    r.source('SMLmetapADJ.R')
    print('called')
    print("reFABSp(c(" + str(pvals)[1:-1] + "))")
    return r("reFABSp(c(" + str(pvals)[1:-1] + "))")

# r = robjects.r
x = reFABScalc([.4,.3,.12])

print(x)