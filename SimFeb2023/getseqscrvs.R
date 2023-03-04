

## ----eval=FALSE---------------------------------------------------------------
## devtools::install_github("https://github.com/simrvprojects/SimRVSequences")


## -----------------------------------------------------------------------------
library(SimRVSequences)
load("Chromwide.Rdata")


## -----------------------------------------------------------------------------
chr8 <- out[[8]]


## -----------------------------------------------------------------------------
library(Matrix) # for sparse matrix functions
keep <- chr8$Mutations$afreq < 0.001 | chr8$Mutations$afreq > 0.999
chr8 <- list(Haplotypes=chr8$Haplotypes[,keep],
              Mutations=chr8$Mutations[keep,])


## -----------------------------------------------------------------------------
# Extract the SNVs from the two TNF genes
posn <- chr8$Mutations$position
TNFmutns <- (posn >= 23020133 & posn <= 23069031) | # TNFRSR10B
            (posn >= 118923557 & posn <= 118951885) # TNFRSR11B
hh <- chr8$Haplotypes[,TNFmutns] #subset haplotypes to area of interest.
mm <- chr8$Mutations[TNFmutns,] #subset mutations to area of interest.
# Create a sampling weight vector for sampling RVs in propn to their
# selection coefficient, stored in the selCoef column.
mm$weight <- abs(mm$selCoef)/sum(abs(mm$selCoef))
# Set a seed for reproducibility and sample RVs with population derived allele count <= 10
# or >= (2*n-10) where n is the number of individuals in the population
# (i.e. twice the number of rows of the dataframe hh)
set.seed(123)
cs <- colSums(hh)
ind <- (cs <= 10 | cs >= (2*nrow(hh) - 10))
ss <- sample((1:nrow(mm))[ind],size=10,prob = mm$weight[ind])
cRVinds <- sort(ss)
cRVschr8 <- mm[cRVinds,]
cat("selected cRVs have total probability:\n")
sum(cRVschr8[,"afreq"])
cat("selected cRVs are:\n")
cRVschr8[,c("position","afreq","selCoef","weight")]
chr8$Mutations$is_CRV <- FALSE # initialize is_CRV col of mutations dataframe
setT <- chr8$Mutations$position %in% cRVschr8$position
chr8$Mutations$is_CRV[setT] <- TRUE # set sampled cRVs to TRUE


## -----------------------------------------------------------------------------
chr8 <- SNVdata(Haplotypes=chr8$Haplotypes,Mutations=chr8$Mutations)
save(chr8,file="chr8.RData")



