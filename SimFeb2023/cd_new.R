cd_new <-  function(peds, subtypes,carrier_probs,
                    tau_increment=0.05,subtype_weights=NULL,useK=FALSE){


  tau_grid <- RVMethods:::make_tauGrid(increment_width = tau_increment, constrained = TRUE)
  study_FamIDs <- unique(peds$FamID)
  famIndex <- c(1:length(study_FamIDs)) # indices of family IDs
  fam_likeGrids <- list() # empty list to hold grids of config probs for each family
  for (i in 1:length(study_FamIDs)){ # loop over families
    ped <- peds[peds$FamID == study_FamIDs[[i]], ]
    # Unlike compute_distributions which has one carrier_prob, we have a vector.
    # Use the first carrier_probs value to compute unconditional familial probs.
    fam_likeGrids[[i]] <- RVMethods:::test_allcombos_withZero(ped, subtypes, tau_grid,
      carrier_probs[1],subtype_weights)
    # Now use our new uncond2cond() function to convert the unconditional
    # familial probs to conditional ones that don't depend on the carrier prob.
    fam_likeGrids[[i]] <- uncond2cond(fam_likeGrids[[i]],ped,carrier_probs[1])
  }

  fam_likeGrids <- lapply(fam_likeGrids, function(x){
    RVMethods:::remove_invalidConfigs(x)
  })

  D1_global_sharing_byBinID <- list() # list for global LR, will be one list item per p_c
  D2_global_sharing_byBinID <- list() # list for global transm,  ditto
  flg <- list()
  for(i in 1:length(carrier_probs)){
    for(j in 1:length(fam_likeGrids)){
      ped <- peds[peds$FamID == study_FamIDs[[j]], ]
      # Need to convert conditional probs in fam_likeGrids to
      # unconditional using current value of carrier_probs before
      # passing to functions that calculate global statistics.
      flg[[j]] <- cond2uncond(fam_likeGrids[[j]],ped,carrier_probs[i])
    }
    D1_global_sharing_byBinID[[i]] <-
      RVMethods:::condition_globalDist_zeroConfig(likeGrids_byFam = flg,
                                                            famID_index = famIndex,
                                                            tau_grid)
    D2_global_sharing_byBinID[[i]] <-
      RVMethods:::approx_globalDist(likeGrids_byFam = flg,
                                               famID_index = famIndex,
                                               tau_grid)

  }
  # local approaches don't depend on p_c, so it doesn't matter
  # what value of the carrier probability was used to get the
  # unconditional familial probs. Can use last val of flg from
  # the for loop above.
  D3_global_sharing_byBinID <-
    RVMethods:::conditioned_semiGlobalDist(likeGrids_byFam = flg,
                                                       famID_index = famIndex,
                                                       tau_grid)

  fam_configs <- list()
  # Find configurations for each family:
  for(i in famIndex){
    fam_configs[[i]] <- 1*fam_likeGrids[[i]]$configs[match(D3_global_sharing_byBinID[, i],
                                                           fam_likeGrids[[i]]$binIDs), ]
    # label columns familyID:subjectID
    colnames(fam_configs[[i]]) <- paste0(study_FamIDs[[i]], ":", colnames(fam_configs[[i]]))
  }
  # fam_configs is now a list of familial configs for each pedigree in the study.
  # Combine them (i.e., cbind() them) into one matrix. The rows of this matrix
  # are the global configs and the columns are the affecteds in the study.
  fam_configs <- do.call(cbind, fam_configs)

  # Global LR: collect LR stats, null probs and calculate p-values.
  # First initialize w/ fam configs, taus and K.
  global_dist <- cbind(fam_configs,
                       D1_global_sharing_byBinID[[1]][, c("tau_A","tau_B")],
                       K=apply(fam_configs, 1, sum))
  for(i in 1:length(carrier_probs)) {
    # add LR and null_configProb cols from ith D1_global_sharing_byBinID matrix
    global_dist <- cbind(global_dist,
                         D1_global_sharing_byBinID[[i]][,c("LR","null_configProb")])
    LR <- paste0("LR",i); nullconf <- paste0("null_configProb",i)
    names(global_dist)[ncol(global_dist) - (1:0)] <- c(LR,nullconf)
    # Calculate p-values and add to global_dist
    if(useK){
      global_dist <- cbind(global_dist,sapply(1:nrow(global_dist), function(x){
        sum(global_dist[global_dist[,LR] >= global_dist[x,LR]
                          & global_dist$K >= global_dist$K[x],nullconf],
            na.rm = TRUE)}))
    } else {
      global_dist <- cbind(global_dist,sapply(1:nrow(global_dist), function(x){
        sum(global_dist[global_dist[,LR] >= global_dist[x,LR],nullconf],
            na.rm = TRUE)}))

    }
    # add index i to LR_pvalue column name
    names(global_dist)[ncol(global_dist)] <- paste0("LR_pvalue",i)
  }

  # Transmission stat: collect LR stats, null probs and calculate p-values.
  # Initialize w/ fam configs, LR (same for all carrier probs), taus and K
  semiglobal_dist <- cbind(fam_configs,
                       D2_global_sharing_byBinID[[1]][, c("LR","tau_A","tau_B")],
                       K=apply(fam_configs, 1, sum))
  for(i in 1:length(carrier_probs)) {
    # add null_configProb col from ith D2_global_sharing_byBinID matrix
    semiglobal_dist <- cbind(semiglobal_dist,
                             D2_global_sharing_byBinID[[i]][,"null_configProb"])
    nullconf <- paste0("null_configProb",i)
    names(semiglobal_dist)[ncol(semiglobal_dist)] <- nullconf
    # Calculate p-values
    if(useK){
      semiglobal_dist$LR_pvalue <- sapply(1:nrow(semiglobal_dist), function(x){
        sum(semiglobal_dist[semiglobal_dist[,"LR"] >= semiglobal_dist[x,"LR"]
                          & semiglobal_dist$K >= semiglobal_dist$K[x],nullconf],
            na.rm = TRUE)})
    } else {
       semiglobal_dist$LR_pvalue <- sapply(1:nrow(semiglobal_dist), function(x){
        sum(semiglobal_dist[semiglobal_dist[,"LR"] >= semiglobal_dist[x,"LR"],nullconf],
            na.rm = TRUE)})

    }
    names(semiglobal_dist)[ncol(semiglobal_dist)] <- paste0("LR_pvalue",i)
  }


  # Local approaches: collect LR stats, null probs and calculate p-values.
  # The following is CN's code
  condsemiglobal_dist <- cbind(fam_configs,
                               D3_global_sharing_byBinID[, -c(1:length(study_FamIDs))])

  condsemiglobal_dist$K = apply(fam_configs, 1, sum)

  #create a list of names of individuals with the genetically compelling subtype
  GCsub_list <- apply(peds[which(peds$available & peds$subtype == subtypes[[1]]), c("FamID", "ID")],
                      1, function(x){paste0(x, collapse = ":")})

  condsemiglobal_dist$K_sub <- sapply(1:nrow(fam_configs), function(x){
    sum(colnames(fam_configs)[fam_configs[x, ] == 1] %in% GCsub_list)
  })

#  print(paste0("Local P-value: ", Sys.time()))
  if(useK){
    #calculate the p-value for the likelihood ratio statistic
    condsemiglobal_dist$LR_pvalue <-
      sapply(1:nrow(condsemiglobal_dist),
      function(x){
      sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$LR >=
                                                condsemiglobal_dist$LR[x]
                                              & condsemiglobal_dist$K >=
                                                condsemiglobal_dist$K[x]
                                              & condsemiglobal_dist$distID ==
                                                condsemiglobal_dist$distID[x]],
          na.rm = TRUE)
    })
  } else {
    condsemiglobal_dist$LR_pvalue <-
      sapply(1:nrow(condsemiglobal_dist),
             function(x){
      sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$LR >=
                                                condsemiglobal_dist$LR[x]
                                              & condsemiglobal_dist$distID ==
                                                condsemiglobal_dist$distID[x]],
          na.rm = TRUE)
    })
  }


#  print(paste0("RVS-Based P-value: ", Sys.time()))
  condsemiglobal_dist$RVS_pvalue <-
    sapply(1:nrow(fam_configs), function(x){
    sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$null_configProb <=
                                              condsemiglobal_dist$null_configProb[x]
                                    & condsemiglobal_dist$K >=
                                      condsemiglobal_dist$K[x]
                                    & condsemiglobal_dist$distID ==
                                      condsemiglobal_dist$distID[x]],
        na.rm = TRUE)
  })

  condsemiglobal_dist$modRVS_pvalue <-
    sapply(1:nrow(fam_configs),
           function(x){
    sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$null_configProb <=
                                              condsemiglobal_dist$null_configProb[x]
                                    & condsemiglobal_dist$K >=
                                      condsemiglobal_dist$K[x]
                                    & condsemiglobal_dist$K_sub >=
                                      condsemiglobal_dist$K_sub[x]
                                    & condsemiglobal_dist$distID ==
                                      condsemiglobal_dist$distID[x]],
        na.rm = TRUE)
  })


  global_dist$binID <-
    apply(global_dist[, 1:ncol(fam_configs)], 1,
          function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  semiglobal_dist$binID <-
    apply(semiglobal_dist[, 1:ncol(fam_configs)], 1,
          function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  condsemiglobal_dist$binID <-
    apply(condsemiglobal_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  global_stats <- paste0("LR",1:length(carrier_probs))
  global_pvals <- paste0("LR_pvalue",1:length(carrier_probs))
    statvals <- cbind(global_dist$binID,  # same ID in all three d.f.s
                    global_dist[,global_stats],
                    semiglobal_dist[,"LR"], # same for all carrier probs
                    condsemiglobal_dist[,"LR"],
                    1/condsemiglobal_dist[,"RVS_pvalue"], # inverse of prob
                    1/condsemiglobal_dist[,"modRVS_pvalue"]) # inverse of prob
  names(statvals) <- c("binID",paste0("globalLR",1:length(carrier_probs)),
                       "globaltransLR","localLR","RVS","modRVS")
  pvals <- cbind(global_dist$binID,  # same ID in all three d.f.s
                 global_dist[,global_pvals],
                 semiglobal_dist[,global_pvals],
                 condsemiglobal_dist[,"LR_pvalue"],
                 condsemiglobal_dist[,"RVS_pvalue"],
                 condsemiglobal_dist[,"modRVS_pvalue"])
  names(pvals) <- c("binID",paste0("globalLR",1:length(carrier_probs)),
                       paste0("globaltrans",1:length(carrier_probs)),
                       "localLR","RVS","modRVS")
  # More convenient for extracting rows to have the lookup tables
  # be matrices rather than data frames.
  output <- list(statvals=as.matrix(statvals), pvals = as.matrix(pvals))

  return(output)
}

# install.packages("devtools")
# library(devtools)
# install_github("https://github.com/simrvprojects/SimRVPedigree")
# install_github("https://github.com/simrvprojects/SimRVSequences")
# install.packages("BiocManager")
# BiocManager::install("RBGL")
# install.packages("gRain")
# install_github("https://github.com/simrvprojects/RVMethods")
# library(RVMethods)
# data(study_pedigrees)
# library(SimRVPedigree)
# # Set up for a function call you can use for testing cd_new()
# peds = study_pedigrees[study_pedigrees$FamID %in% c(58, 304), ]
# subtypes = c("HL", "NHL")
# carrier_prob = 0.00032 # True carrier prob
# carrier_probs = carrier_prob*c(1/10,1/2,1,2,10)
# # Call cd_new()
# cd_new(peds,subtypes,carrier_probs) #leave remaining arguments at their defaults

uncond2cond <- function(fam_likeGrid,ped,carrier_prob) {
  like_mat <- fam_likeGrid$like_mat
  num_found = length(unique(ped$ID[which(is.na(ped$dadID) & is.na(ped$momID))]))
  fint_prob = num_found*carrier_prob #founder intro prob, using first carrier prob
  like_mat[,1:(ncol(like_mat)-1)] <-
      like_mat[,1:(ncol(like_mat)-1)] / fint_prob #non-zero configs
  like_mat[,ncol(like_mat)] <-
      (like_mat[,ncol(like_mat)] - (1-fint_prob))/fint_prob # zero config
  fam_likeGrid$like_mat <- like_mat
  return(fam_likeGrid)
}
# Also do the reverse function
cond2uncond <- function(fam_likeGrid,ped,carrier_prob) {
  like_mat <- fam_likeGrid$like_mat
  num_found = length(unique(ped$ID[which(is.na(ped$dadID) & is.na(ped$momID))]))
  fint_prob = num_found*carrier_prob #founder intro prob, using first carrier prob
  like_mat[,1:(ncol(like_mat)-1)] <-
      fint_prob*like_mat[,1:(ncol(like_mat)-1)] # non-zero configs
  like_mat[,ncol(like_mat)] <-
      fint_prob*like_mat[,ncol(like_mat)] + (1-fint_prob) # zero config
  fam_likeGrid$like_mat <- like_mat
  return(fam_likeGrid)
}
