This folder contains scripts to resurrect the simulation study in Christina's thesis.

1. To get started on simulating 150 pedigrees, see:

* **simRVped.Rmd**: An Rmarkdown file describing how to simulate pedigrees in the R script `simrvped.R`. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **simrvped.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **simrvped.sh**: The SLURM script to run an array job on the cluster.  To run the script, you can type "SBATCH simrvped.sh" from the command line on the cluster (appendix A of the group's Workflow document may be helpful to fill in more details).

Once we get the 150 pedigrees, we screen them according to the criteria in Christina's thesis. She only used 55 pedigrees. The screening:
* checks the 150 simulated pedigrees against the selection criteria in the simulation study of Christina's thesis and 
* obtains 55 of the eligible pedigrees to work with. 

The pedigree-screening files are:
* **checkpeds.Rmd**: An Rmarkdown file describing how the R script checkpeds.R checks the simulated pedigrees and obtains the 55 we will work with. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **checkpeds.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **checkpeds.sh**: The SLURM script to run on the cluster (not strictly necessary).  To run the script, you can type "SBATCH checkpeds.sh" from the command line on the cluster.

Once we get the 55 eligible pedigrees to work with, the next step will be to "seed" their founders with exome sequences from Nirodha's simulated American admixed population and specify the causal RVs (cRVs) for the "conditional gene drop" in SimRVPedigree.
