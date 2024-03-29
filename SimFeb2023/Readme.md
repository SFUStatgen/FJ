This folder contains scripts to resurrect the simulation study in Christina's thesis. An overview of the simulation workflow is provided in the file **simoverview.pdf**, which is generated by the Rmarkdown file **simoverview.Rmd**. Tasks in the simulation workflow can have up to three associated scripts: an Rmarkdown file, an R script and a SLURM script. The Rmarkdown file is central because it generates a non-working draft of the R script as well as PDF documentation for the task. The SLURM script, when necessary, is provided separately and runs the R script on the cluster. The idea is to run the Rmarkdown file on your PC to generate a draft version of the R script and then edit the R script so that it works (usually by removing the first and last commands in the file). The PDF documentation for the task will be also rendered by running the Rmarkdown file on your PC. Afterwards, you can port all the scripts from your PC to the cluster (as well as the PDF documentation if you want). 

1. After reading through **simoverview.pdf**, you can get started on simulating 150 pedigrees with:

* **simRVped.Rmd**: An Rmarkdown file describing how to simulate pedigrees in the R script `simrvped.R`. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **simrvped.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **simrvped.sh**: The SLURM script to run an array job on the cluster.  To run the script, you can type "SBATCH simrvped.sh" from the command line on the cluster (appendix A of the group's Workflow document may be helpful to fill in more details).

Once we get the 150 pedigrees, we screen them according to the criteria in Christina's thesis. She used only 55 pedigrees in her pool of pedigrees to sample from to make a "study" (comprised of three pedigrees). The screening:
* checks the 150 simulated pedigrees against the selection criteria in the simulation study of Christina's thesis and 
* obtains 55 of the eligible pedigrees to work with. 

The pedigree-screening files are:
* **checkpeds.Rmd**: An Rmarkdown file describing how the R script **checkpeds.R** checks the simulated pedigrees and obtains the 55 we will work with. This Rmarkdown file assumes that you have installed the SimRVPedigree package.
* **checkpeds.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the SimRVPedigree package in your account on the cluster.
* **checkpeds.sh**: The SLURM script to run on the cluster (not strictly necessary).  To run the script, you can type "SBATCH checkpeds.sh" from the command line on the cluster.

2. Once we get the 55 eligible pedigrees to work with, the next step is to get the pool of chromosome 8 exome sequences from which to draw pedigree founder sequences from Nirodha's simulated American admixed population and specify the pool of causal RVs (cRVs) in the population. The relevant files for this step are:
* **getseqscrvs.Rmd**: An Rmarkdown file describing how the R script **getseqscrvs.R** gets the chromosome 8 sequences and cRVs in the population. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository (https://zenodo.org/record/6369360#.ZEfwRezMI6E) and placed in the same directory.
* **getseqscrvs.R**: The associated R script called by the SLURM script on the Compute Canada cluster. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository and placed in the same directory.
* **getseqscrvs.sh**: The SLURM script to run on the cluster (not strictly necessary).  To run the script, you can type "SBATCH getseqscrvs.sh" from the command line on the cluster. Assumes that the R data file chromwide.Rdata has been downloaded from Nirodha's Zenodo repository and placed in the same directory.

3. Once we have the chromosome 8 founder sequences for the 55 pedigrees we can start running simulations to get the power and ranking results under the alternative hypothesis. The files for running simulations under the alternative hypothesis are:
* **simalt.Rmd**: An Rmarkdown file describing how the R script **simalt.R** does the simulations under the alternative hypothesis. This Rmarkdown file assumes that you have installed the RVMethods package and that the 55 eligible pedigrees and the pool of chromosome 8 exome sequences are on the Compute Canada cluster in the `/project/def-jgraham/FJdata` directory. 
* **simalt.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the RVMethods package in your account on the cluster and that the 55 eligible pedigrees and the pool of chromosome 8 exome sequences are on the Compute Canada cluster in the `/project/def-jgraham/FJdata` directory. 
    * The R script calls a function called `cd_new()` to get lookup tables of p-values and test statistics for the five methods we consider. `cd_new()` is based on the function `compute_distributions()` from RVMethods, with modifications to make it more efficient for simulations. The function is documented in the Rmarkdown file **cd_new.Rmd**, with associated R script **cd_new.R**.
* **simalt.sh**: The SLURM script to run on the cluster.  To run the script, you can type "SBATCH simalt.sh" from the command line on the cluster.

4. When the simulations under the alternative hypothesis finish on the cluster, copy the output files `pvalres`**i**.`csv` and `rankres`**i**.`csv` for **i**`=1..200`, back to your PC and prepare the summaries for the manuscript by knitting the Rmarkdown file **simaltSummary.Rmd** on your PC.

5. The files for running simulations under the null hypothesis are:
* **simnull.Rmd**: An Rmarkdown file describing how the R script **simnull.R** does the simulations. This Rmarkdown file assumes that you have installed the RVMethods package and that the 55 eligible pedigrees and the pool of chromosome 8 exome sequences are on the Compute Canada cluster in the `/project/def-jgraham/FJdata` directory. 
* **simnull.R**: The associated R script called by the SLURM script on the Compute Canada cluster. This R script assumes that you have installed the RVMethods package in your account on the cluster and that the 55 eligible pedigrees and the pool of chromosome 8 exome sequences are on the Compute Canada cluster in the `/project/def-jgraham/FJdata` directory. 
    * The R script calls a function called `cd_new()` to get lookup tables of p-values and test statistics for the five methods we consider. `cd_new()` is based on the function `compute_distributions()` from RVMethods, with modifications to make it more efficient for simulations. The function is documented in the Rmarkdown file **cd_new.Rmd**, with associated R sript **cd_new.R**.
* **simnull.sh**: The SLURM script to run on the cluster.  To run the script, you can type "SBATCH simnull.sh" from the command line on the cluster.

6. When the simulations under the null hypothesis finish on the cluster, copy the output files `pvalnullires`**i**.`csv` for **i**`=1..200` back to your PC and prepare the summaries for the manuscript by knitting the Rmarkdown file **simnullSummary.Rmd** on your PC.
