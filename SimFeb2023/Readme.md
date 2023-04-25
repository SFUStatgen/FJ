This folder contains scripts to resurrect the simulation study in Christina's thesis. Each task has three scripts associated with it: an Rmarkdown file, an R script and a SLURM script. The Rmarkdown file documents the commands in the R script and the SLURM script runs the R script on the cluster. The idea is to run the Rmarkdown file on your PC to generate a draft version of the R script and then edit the R script so that it works (usually by removing the first and last commands in the file). Afterwards, you can port all the scripts, including the PDF documentation rendered by running the Rmarkdown file on your PC, to the cluster. 

1. To get started on simulating 150 pedigrees, see:

* **simRVped.Rmd**: An Rmarkdown file describing how to simulate pedigrees in the R script `simrvped.R`. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **simrvped.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **simrvped.sh**: The SLURM script to run an array job on the cluster.  To run the script, you can type "SBATCH simrvped.sh" from the command line on the cluster (appendix A of the group's Workflow document may be helpful to fill in more details).

Once we get the 150 pedigrees, we screen them according to the criteria in Christina's thesis. She only used 55 pedigrees. The screening:
* checks the 150 simulated pedigrees against the selection criteria in the simulation study of Christina's thesis and 
* obtains 55 of the eligible pedigrees to work with. 

The pedigree-screening files are:
* **checkpeds.Rmd**: An Rmarkdown file describing how the R script **checkpeds.R** checks the simulated pedigrees and obtains the 55 we will work with. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **checkpeds.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **checkpeds.sh**: The SLURM script to run on the cluster (not strictly necessary).  To run the script, you can type "SBATCH checkpeds.sh" from the command line on the cluster.

2. Once we get the 55 eligible pedigrees to work with, the next step is to get the pool of exome sequences from which to draw pedgree founder sequences from Nirodha's simulated American admixed population and specify the pool of causal RVs (cRVs) in the population. The relevant files for this step are:
* **getseqscrvs.Rmd**: An Rmarkdown file describing how the R script **getseqscrvs.R** gets the sequences and cRVs in the population. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository (https://zenodo.org/record/6369360#.ZEfwRezMI6E) and placed in the same directory.
* **getseqscrvs.R**: The associated R script called by the SLURM script on the Compute Canada cluster. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository and placed in the same directory.
* **getseqscrvs.sh**: The SLURM script to run on the cluster (not strictly necessary).  To run the script, you can type "SBATCH getseqscrvs.sh" from the command line on the cluster. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository and placed in the same directory.

3. Once we get the 55 pedigrees and the founder sequences we can start doing our simulations. The files for doing simulations under the alternative hypothesis are:
* **simalt.Rmd**: An Rmarkdown file describing how the R script **simalt.R** does the simulations. This Rmarkdown file assumes that you have installed the RVMethods package. 
* **simalt.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the RVMethods package in your account on the cluster. The R script calls a function called `cd_new()` to get lookup tables of p-values and test statistics for the five methods we consider. `cd_new()` is based on the function `compute_distributions()` from RVMethods, with modifications to make it more efficient for simulations. The function is documented in the Rmarkdown file **cd_new.Rmd**, with associated R sript **cd_new.R**.
* **simalt.sh**: The SLURM script to run on the cluster.  To run the script, you can type "SBATCH simalt.sh" from the command line on the cluster.
