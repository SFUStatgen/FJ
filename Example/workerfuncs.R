## -----functions written for exdatan.Rmd and documented in the-------------
## -----Appendix of that .Rmd file------------------------------------------
remove_invalid_configs <- function(a_seqs,lookupTabs) {
  # configs are identified in the lookup table by their base-10
  # representation, or "binID"
  valid_binIDs <- lookupTabs$statvals[,"binID"]
  # initialize an empty vector to hold info on which RV configs are valid
  include <- rep(NA,nrow(a_seqs$SNV_map))
  # loop over RVs
  for(i in 1:length(include)) {
    RV_binID <- find_binID(a_seqs$SNV_map$marker[i],a_seqs)
    include[i] <- (RV_binID %in% valid_binIDs)
  }
  # Use the include vector to subset the input a_seqs and return
  out <- list(ped_files=a_seqs$ped_files,
              ped_haplos=a_seqs$ped_haplos[,include],
              haplo_map=a_seqs$haplo_map,
              SNV_map=a_seqs$SNV_map[include,],
              ped_genos=a_seqs$ped_genos[,include])
  return(out)
}
# Modified version of Christina's summary.famStudy() function. The modification
# is to split allele counts by subtype. Very quick-and-dirty, and assumes
# three families in the study with 12 affecteds in total.
mysummary.famStudy <- function (object, ...)
{
  Fids <- sort(unique(object$ped_files$FamID))
  # Total allele count in affecteds
  aff_allele_counts <- lapply(Fids, function(x) {
    SimRVSequences:::affected_allele_count(
      ped_haps = object$ped_haplos[object$haplo_map$FamID == x, ],
      hap_map = object$haplo_map[object$haplo_map$FamID == x, ],
      ped_file = object$ped_files[object$ped_files$FamID == x, ])
  })
  tot_allele_count <- do.call(rbind, aff_allele_counts)
  totals <- colSums(tot_allele_count)
  # Alelle counts in those affected with NHL
  aff_allele_counts <- lapply(Fids, function(x) {
    affected_allele_count_NHL(
      ped_haps = object$ped_haplos[object$haplo_map$FamID == x, ],
      hap_map = object$haplo_map[object$haplo_map$FamID == x, ],
      ped_file = object$ped_files[object$ped_files$FamID == x, ])
  })
  NHL_allele_count <- do.call(rbind, aff_allele_counts)
  # Alelle counts in those affected with HL
  aff_allele_counts <- lapply(Fids, function(x) {
    affected_allele_count_HL(
      ped_haps = object$ped_haplos[object$haplo_map$FamID == x, ],
      hap_map = object$haplo_map[object$haplo_map$FamID == x, ],
      ped_file = object$ped_files[object$ped_files$FamID == x, ])
  })
  HL_allele_count <- do.call(rbind, aff_allele_counts)
  fam_allele_count <- cbind(
    NHL_allele_count[ 1,],HL_allele_count[1,],tot_allele_count[1,],
    NHL_allele_count[2,],HL_allele_count[2,],tot_allele_count[2,],
    NHL_allele_count[3,],HL_allele_count[3,],tot_allele_count[3,]
  )
  colnames(fam_allele_count) <- c("NHL1","HL1","tot1","NHL2","HL2","tot2","NHL3","HL3","tot3")
  rownames(fam_allele_count) <- object$SNV_map$marker
  out <- data.frame(
    chrom = object$SNV_map$chrom,
    position = object$SNV_map$position,
    fam_allele_count)
  return(out[totals>0 & totals < 24,]) # 24 is count for all carriers of 1 allele
}

affected_allele_count_NHL <- function (ped_haps, hap_map, ped_file)
{
  aff_IDs <- ped_file$ID[ped_file$affected & ped_file$subtype=="NHL"]
  aff_rows <- which(hap_map$ID %in% aff_IDs)
  total_count <- colSums(ped_haps[aff_rows, ])
  return(total_count)
}
affected_allele_count_HL <- function (ped_haps, hap_map, ped_file)
{
  aff_IDs <- ped_file$ID[ped_file$affected & ped_file$subtype=="HL"]
  aff_rows <- which(hap_map$ID %in% aff_IDs)
  total_count <- colSums(ped_haps[aff_rows, ])
  return(total_count)
}

## --------------------------------------------------------------------------
# The functions from here on are taken from R scripts in the simulation
# workflow and are documented in the corresponding .Rmd files.
## --------------------------------------------------------------------------

## -------functions from checkpeds.R - see checkpeds.Rmd for documentation----
plot_affped <- function(ped,cex=.4,ref_year=2018){
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
  ped <- ped[ped$ID %in% keepIDs,]
  # Use Christina's plot function on the reduced pedigree
  SimRVPedigree:::plot.ped(ped,cex=cex,ref_year)
  invisible(ped)
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

## -----functions from simalt.R - see simalt.Rmd for documentation---------

# read_studypeds() takes a vector of pedigree IDs and reads the
# pedigrees from their plain-text files
read_studypeds <- function(studypedIDs,infileDir){
  for(i in 1:length(studypedIDs)){
    pedfile <- paste0(infileDir,"/ascertained_ped",studypedIDs[i],".txt")
    # Read in pedfile and set the resulting object to be of class "ped" and "data.frame".
    pp <- read.table(pedfile);class(pp) <- c("ped","data.frame")
    # Set FamID column to the ID of our sampled pedigree.
    pp[,"FamID"] <- studypedIDs[i]
    if(i==1) { # Initialize the data frame containing the pedigrees.
      s_peds <- pp
    } else{ # Add new pedigree to the data frame of previous ones.
      s_peds <- rbind(s_peds,pp)
    }
  }
  return(s_peds) #Return the pedigrees.
}
# filter2aff() takes a "famStudy" object returned by sim_RVstudy()
# and reduces its sequence data to rare variants that are seen in
# the affected individuals. It also finds the genotypes for each
# variant/individual. We'll need these genotypes later in the
# ranking application, when we find the RV configurations and their
# corresponding test stats for each variant in the sequence data.
filter2aff <- function(seqs){
  # seqs is a famStudy object, which is a list having elements
  # - ped_haplos (sparse matrix of sequence data),
  # - haplo_map (maps the sequences to individuals), and
  # - SNV_map (support info on each SNV).
  # First reduce to sequences in affected pedigree members with DNA
  ped_haplos <- seqs$ped_haplos[seqs$haplo_map$affected,]
  haplo_map <- seqs$haplo_map[seqs$haplo_map$affected,]
  # Next filter variants to those that appear in the affected individuals with DNA.
  # The variants appear as columns of ped_haplos, and rows of SNV_map.
  cc <- colSums(ped_haplos)
  ped_haplos <- ped_haplos[,cc>0]
  SNV_map <- seqs$SNV_map[cc>0,]
  # Lastly, pair sequences from individuals into multilocus genotypes
  odd_inds <- seq(from=1,to=nrow(ped_haplos)-1,by=2)
  even_inds <- seq(from=2,to=nrow(ped_haplos),by=2)
  ped_genos <- ped_haplos[odd_inds,] + ped_haplos[even_inds,]
  return(list(ped_files=seqs$ped_files,ped_haplos=ped_haplos,
              haplo_map=haplo_map,
              SNV_map=SNV_map, ped_genos=ped_genos))
}
# find_binID takes the name of an RV and the genotypes of the affected
# individual returned by filter2aff() and returns the compact base-10
# representation of the configuration (what Christina calls the "binID").
find_binID <- function(RV,a_seqs) {
  config_vec <- a_seqs$ped_genos[,a_seqs$SNV_map$marker == RV] # vec of 0's and 1's
  return(config2binID(config_vec))
}
config2binID <- function(config_vec) {
  # Collapse the vector to a string of 0's and 1's and then use
  # strtoi to convert this base-2 number to an integer.
  return(base::strtoi(paste0(config_vec, collapse = ""), base = 2))
}
# get_seqstats takes the observed sequence data on affected individuals
# returned by filter2aff() and the lookup table of statistic values
# for all global configurations and returns statistic values for each
# configuration observed in the sequence data.
get_seqstats <- function(ped_genos,statvals){
  # - ped_genos is the genotypes of the affected individuals at the
  # variants that are observed in the affected individuals.
  # - statvals is the lookup table of statistic values
  binID <- statvals[,"binID"] # binIDs of rows of statvals
  allstats <- matrix(NA,nrow=ncol(ped_genos),ncol=ncol(statvals))
  colnames(allstats) <- colnames(statvals)
  for(i in 1:ncol(ped_genos)) {
    binID.obs <- config2binID(ped_genos[,i]) # *worker*
    if(binID.obs %in% binID) { # then this config is in the lookup table
      allstats[i,1] <- binID.obs
      allstats[i,] <- statvals[binID==binID.obs,]
    }
  } # End loop over the SNVs.
  foundconfig <- !is.na(allstats[,1])
  allstats <- allstats[foundconfig,]
  return(allstats)
}
# get_pedcRVs finds the cRVs of each pedigree in the study
get_pedcRVs <- function(a_seqs){
  # Start by tabulating unique combinations of the family
  # ID and cRV ID. Then extract the cRV IDs.
  FamIDcRV <- unique(a_seqs$haplo_map[,c("FamID","FamCRV")])
  return(FamIDcRV$FamCRV)
}
get_pvals <- function(cRVs,lookupTabs,a_seqs){
  # Create a matrix to hold the p-values, with rows for cRVs and columns
  # for the different methods. In this context, the global LR and global
  # transmission tests with different carrier probs are different "methods".
  # The number of methods is the number of columns of lookupTabs$pvals
  # minus 1 (since the first column of lookupTabs$pvals is binID).
  pvals <- matrix(NA,nrow=3,ncol=(ncol(lookupTabs$pvals)-1))
  for(i in 1:length(cRVs)) {
    # Find binID of the current cRV's configuration.
    binID <- find_binID(cRVs[i],a_seqs) # *worker*
    # Extract pvals for this configuration from the lookup table but
    # remember that the first column of the lookup table is the binID.
    pvals[i,] <- lookupTabs$pvals[lookupTabs$pvals[,"binID"]==binID,-1]
  }
  return(pvals)
}
get_ranks <- function(cRVs,lookupTabs,a_seqs){
  # Create a matrix to hold the ranks, with rows for cRVs and columns for
  # the different methods. In this context, the global LR approach with
  # different values of the carrier probability is taken to be different
  # "methods". The # of methods is the # of columns of lookupTabs$statvals
  # minus 1 (since the 1st column of lookupTabs$statvals is binID).
  ranks <- matrix(NA,nrow=max(3,length(cRVs)),ncol=(ncol(lookupTabs$statvals)-1))
  for(i in 1:length(cRVs)) {
    # Find base-10 representation, or "binID", of the configuration.
    binID <- find_binID(cRVs[i],a_seqs)
    # Calculate ranks over the values of the statistics corresponding to
    #the observed configurations in the sequence data on affected individuals.
    obs_stats <- get_seqstats(a_seqs$ped_genos,lookupTabs$statvals) # *worker*
    for(j in 1:ncol(ranks)){
      # Rank current cRV on jth statistic (in (j+1)st col of obs_stats)
      curr_stat <- unique(obs_stats[obs_stats[,"binID"]==binID,j+1])
      ranks[i,j]<-mean(obs_stats[,j+1] >= curr_stat)
    }
  }
  return(ranks)
}

