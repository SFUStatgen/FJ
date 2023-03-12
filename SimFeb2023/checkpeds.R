## -----------------------------------------------------------------------------
plot_affped <- function(ped,cex=.4){
  # Get IDs of affecteds in the pedigrees
  aff <- ped$affected # may include missing values
  aff[is.na(aff)] <- FALSE
  affIDs <- ped$ID[aff]
  # Initialize a vector of pedigree member IDs we'd like to keep
  # because they are an affected member or their ancestor
  keepIDs <- NULL
  # Loop over affecteds and add their ID and their ancestors' IDs to keepIDs
  for(i in 1:length(affIDs)){
    keepIDs <- c(keepIDs,ancestorIDs(affIDs[i],ped))
  }
  # remove duplicates
  keepIDs <- unique(keepIDs)
  # Use Christina's plot function on the reduced pedigree
  SimRVPedigree:::plot.ped(ped[ped$ID %in% keepIDs,],cex=cex)
}
ancestorIDs <- function(ID,ped){
  # Start with input (child) ID and then call self on mom and dad.
  # We will stop the recursion when we reach a founder, which has missing
  # mom and dad IDs
  keepIDs <- ID
  momID <- ped[ped$ID==ID,"momID"]
  if(!is.na(momID)) keepIDs <- c(keepIDs,ancestorIDs(momID,ped))
  dadID <- ped[ped$ID==ID,"dadID"]
  if(!is.na(dadID)) keepIDs <- c(keepIDs,ancestorIDs(dadID,ped))
  return(keepIDs)
}


## -----------------------------------------------------------------------------
set.seed(42)
library(SimRVPedigree)
# simulation study parameters
npeds <- 150 # number of simulated pedigrees
npool <- 55 # size of pool of pedigrees we want for our simulation study
# vectors to hold indices of excluded pedigrees
morethan6 <- NULL # more than 6 disease-affected members
noHLcarr <- NULL # No HL carriers of cRV in affecteds
noNHL <- NULL # No NHL among affecteds
for(i in 1:npeds){
  cat("\n----------------Pedigree",i,"-------------------\n")
  pfile <-  paste0("Outputfiles/ascertained_ped",i,".txt")
  if(!file.exists(pfile)) {
    cat("file",pfile,"does not exist!\n")
  } else { # read the pedigree into R and check it
    ped <-  read.table(pfile)
    # ped is read in as a data frame. It also needs to have class "ped"
    # to be plotted.
    class(ped) <-c("ped","data.frame")
    plot_affped(ped,cex=.4)
    # Now start checking the pedigree. First check for missing affection
    # status but non-missing subtype.
    if(any(is.na(ped$affected) & !is.na(ped$subtype))) {
      cat("***ped members with missing affected and non-missing subtype -- \n")
      cat("   setting such subtypes to missing and re-writing ped to file \n")
      is.na(ped$subtype) <- is.na(ped$affected)
      write.table(ped,file=pfile)
    }
    # Check for >6 affecteds.
    pind <- ped$affected
    pind[is.na(pind)] <- FALSE # set missing affection status to FALSE
    if(sum(pind)>6) {
      cat("***More than 6 disease-affecteds***\n")
      morethan6 <- c(morethan6,i)
    }
    # The last two checks pertain to the affecteds, so subset the
    # pedigree to just the affecteds.
    aped <- ped[pind,]
    # Check for no HL carriers of a disease allele
    aped$disease.alleles = aped$DA1+aped$DA2
    if(sum(aped$subtype=="HL" & aped$disease.alleles>0)==0){
      cat("***No HL carriers ***")
      noHLcarr <- c(noHLcarr,i)
    }
    # Lastly, check for no NHLs in the pedigree
    if(!any(aped$subtype=="NHL")) {
      cat("***No NHL subtype***\n")
      noNHL <- c(noNHL,i)
    }
  }
}
# Print out the indices of the pedigrees to be excluded
cat("more than 6 disease-affecteds:\n")
print(morethan6)
cat("disease alleles in NHL, but none in HL:\n")
print(noHLcarr)
cat("no NHL subtype among affecteds:\n")
print(noNHL)
# Keep the pedigrees not in one of the exclusion vectors.
keep <- setdiff(1:npeds,c(morethan6,noHLcarr,noNHL))
cat("keep peds:\n")
print(keep)


## -----------------------------------------------------------------------------
pedpool <- keep[1:npool]
write(pedpool,file="Outputfiles/pedpool.txt")
