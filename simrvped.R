
## ------------------------------------------------------------------------------------------
library(SimRVPedigree)
data(SubtypeHazards)
head(SubtypeHazards)

my_hazards <- hazard(SubtypeHazards,
                     subtype_ID = c("HL", "NHL"))


## ------------------------------------------------------------------------------------------
dID = Sys.getenv("SLURM_ARRAY_TASK_ID")
seed = as.numeric(dID)
# Set a seed value to assure the reproducibility.
if(!is.na(seed)) {
  set.seed(seed)
} else {
  warning("No task ID, setting seed to 1")
  set.seed(1)
}


## ------------------------------------------------------------------------------------------
generatePeds = function(dataID){
  # Simulate pedigree ascertained for at least two individuals
  # affected by either Hodgkin's lymphoma or non-Hodgkin's lymphoma.
  out <- sim_RVped(hazard_rates = my_hazards,
                      GRR = c(35, 1),
                      RVfounder = TRUE,
                      FamID = 1,
                      founder_byears = c(1825, 1850),
                      ascertain_span = c(2000, 2010),
                      num_affected = 4,
                      stop_year = 2018,
                      carrier_prob = 0.00032, # 2x cum prob of cRVs
                      recall_probs = c(1, 1, 1, .75,.125,.125, 0),
                      first_diagnosis = 1980,
                      sub_criteria = list("HL",1)) # ascertain only if at least one HL

   write.table(out$full_ped, file = paste0("Outputfiles/full_ped",dataID,".txt"))
   write.table(out$ascertained_ped, file = paste0("Outputfiles/ascertained_ped",dataID,".txt"))
}
# Run the function.
generatePeds(dID)

