## -------------------------------------------------------------------------------------
infileDir <- "/project/def-jgraham/FJdata"
outfileDir <- "/project/def-jgraham/FJdata"


## -------------------------------------------------------------------------------------
dID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(dID)) {
  dID=1
  warning("No task ID, setting task ID to 1")
}
seed <- dID
set.seed(seed)


## -------------------------------------------------------------------------------------
library(SimRVSequences)
library(RVMethods)
library(Matrix)
# Load the cd_new() function from its R source file
source("cd_new.R")


## -------------------------------------------------------------------------------------
pedpool <- scan(paste0(infileDir,"/pedpool/pedpool.txt"))
load(paste0(infileDir,"/chr8.RData"))


## -------------------------------------------------------------------------------------
oldop <- options(scipen=999)


## -------------------------------------------------------------------------------------
N <- 10  # number of simulation reps per batch
studysize <- 3 # number  of pedigrees in a study


## -------------------------------------------------------------------------------------
true_carrier_prob <- 0.00032
carrier_probs <- true_carrier_prob*c(1/10,1/2,1,2,10)


## -------------------------------------------------------------------------------------
pfile <- paste0(outfileDir,"/pvalnullres",dID,".csv")
repcols <-c("rep","studyped_1","studyped_2","studyped_3","cRV_1","cRV_2","cRV_3")
pvalcols <- c(paste0("globalLR",1:length(carrier_probs),"_1"),
               paste0("globaltrans",1:length(carrier_probs),"_1"),
               "localLR_1","RVS_1","modRVS_1",
               paste0("globalLR",1:length(carrier_probs),"_2"),
               paste0("globaltrans",1:length(carrier_probs),"_2"),
               "localLR_2","RVS_2","modRVS_2",
               paste0("globalLR",1:length(carrier_probs),"_3"),
               paste0("globaltrans",1:length(carrier_probs),"_3"),
               "localLR_3","RVS_3","modRVS_3")
cat(paste0(c(repcols,pvalcols),collapse=","),"\n",file=pfile) # write header to file



## -------------------------------------------------------------------------------------
simnullreps <- function(N,studysize,pfile,infileDir){
for(simrep in 1:N) {
  cat("simrep",simrep,"\n")
  #------------------------------------------------------------
  # a. Sample the IDs of 3 peds from pool of good pedigrees
  studypedIDs <- sample(pedpool,size=studysize)
  # Read the three pedigrees into a single data frame
  s_peds <- read_studypeds(studypedIDs,infileDir) # *worker*
  # Remove the DA1 and DA2 columns of the pedigree data structure so that
  # sim_RVstudy() will do unconditional gene drop.
  s_peds$DA1 <- s_peds$DA2 <- NULL
  #------------------------------------------------------------
  # b. Generate sequence data with sim_RVstudy() until we get one with a
  # candidate cRV that is present in the affected individuals.
  foundPoly <- FALSE; genedropcounter<- 0
  while(!foundPoly) {
    genedropcounter <- genedropcounter+1
    s_seqs <- sim_RVstudy(s_peds,chr8) # *Christina's*
    # see if any cRVs are polymorphic in the affecteds
    if(sum(s_seqs$ped_haplos[s_seqs$haplo_map$affected,s_seqs$SNV_map$is_CRV]) > 0){
      foundPoly <- TRUE
    }
  }
  # Report how many gene drops we had to do
  cat("Did",genedropcounter,"genedrops to get a polymorphic cRV in affecteds \n")
  # Post-simulation filtering to affecteds and SVNs that
  # are polymorphic in the affecteds
  a_seqs <- filter2aff(s_seqs) # *worker*
  numRVs <- nrow(a_seqs$SNV_map)
  #------------------------------------------------------------
  # c. call cd_new() to get lookup tables of statistics and p-values for
  # each possible global configurations of affecteds in the study.
  lookupTabs = cd_new(peds = s_peds, subtypes = c("HL", "NHL"),
                      carrier_probs = carrier_probs)
  # The output lookupTabs is a list with elements statvals and pvals
  # d. For each familial cRV, find its global configuration. Then find
  # pvalues and write these to an output file. In the code below we assume
  # that there are at most three polymorphic markers from our list of
  # candidate cRVs in a given study. (Turns out there was never more than two.)
  pvals <- matrix(NA,nrow=3,ncol=3+2*length(carrier_probs))
  polycRVs <- a_seqs$SNV_map$marker[a_seqs$SNV_map$is_CRV]
  for(i in 1:length(polycRVs)) {
    # find base-10 representation of config (what Christina calls "binID")
    binID <- find_binID(polycRVs[i],a_seqs) # *worker*
    # extract pvals for this config from the lookup tables in lookupTabs
    # jth pval is in column j+1 of pvals lookup table (binID is in 1st col)
    pvals[i,] <- lookupTabs$pvals[lookupTabs$pvals[,"binID"]==binID,-1]
  }
  # Write the pval and rank results to their files. Also write info
  # on the rep number, pedigree IDs and cRVs in each pedigree
  polycRVs <- polycRVs[1:3] # force vector of length 3, with NAs if len < 3
  repinfo <- c(simrep,studypedIDs,polycRVs)
  pvalinfo <- c(pvals[1,],pvals[2,],pvals[3,])
  cat(paste0(paste0(c(repinfo,pvalinfo),collapse=","),"\n"),
      file=pfile,append=TRUE)
} # end for loop over reps
  return(NULL) # output is written to files, so nothing to return
} # end simnullreps()


## -------------------------------------------------------------------------------------
# read_studypeds() takes a vector if pedigree IDs and reads the
# pedigrees from their plain-text files
read_studypeds <- function(studypedIDs,infileDir){
  for(i in 1:length(studypedIDs)){
    pedfile <- paste0(infileDir,"/ascertained_ped",studypedIDs[i],".txt")
     # read in pedfile and set object to be of class "ped" and "data.frame"
    pp <- read.table(pedfile);class(pp) <- c("ped","data.frame")
    # Set FamID column to the ID of our sampled ped
    pp[,"FamID"] <- studypedIDs[i]
    if(i==1) { # initialize output pedigreees data frame
      s_peds <- pp
    } else{ # add new ped to the data frame of previous ones
      s_peds <- rbind(s_peds,pp)
    }
  }
  return(s_peds)
}
# filter2aff() takes a "famStudy" object output by sim_RVstudy() and
# reduces its sequence data to variants that are polymorphic in the
# affected individuals.
filter2aff <- function(seqs){
  # seqs is a list with elements
  # - ped_haplos (sparse matrix of sequence data),
  # - haplo_map (maps the sequences to individuals), and
  # - SNV_map (support info on each SNV).
  # First reduce to sequences in affecteds
  ped_haplos <- seqs$ped_haplos[seqs$haplo_map$affected,]
  haplo_map <- seqs$haplo_map[seqs$haplo_map$affected,]
  # Next filter variants to those that appear in affected individuals.
  # The variants appear as columns of ped_haplos, and rows of SNV_map
  cc <- colSums(ped_haplos)
  ped_haplos <- ped_haplos[,cc>0]
  SNV_map <- seqs$SNV_map[cc>0,]
  # Lastly, pair sequences from individuals into multilocus genotypes
  odd_inds <- seq(from=1,to=nrow(ped_haplos)-1,by=2)
  even_inds <- seq(from=2,to=nrow(ped_haplos),by=2)
  ped_genos <- ped_haplos[odd_inds,] + ped_haplos[even_inds,]
  return(list(ped_haplos=ped_haplos, haplo_map=haplo_map,
              SNV_map=SNV_map, ped_genos=ped_genos))
}
# find_binID takes the name of an RV and the reduced sequence data on
# affecteds output by filter2aff() and returns the base-10 representation
# of the configuration (what Christina calls the "binID").
find_binID <- function(RV,a_seqs) {
  config_vec <- a_seqs$ped_genos[,a_seqs$SNV_map$marker == RV] # vec of 0's and 1's
  return(config2binID(config_vec))
}
config2binID <- function(config_vec) {
  # collapse the vector to a string of 0's and 1's and then use
  # strtoi to convert this base-2 number to an integer.
  return(base::strtoi(paste0(config_vec, collapse = ""), base = 2))
}


## -------------------------------------------------------------------------------------
simnullreps(N,studysize,pfile,infileDir)

