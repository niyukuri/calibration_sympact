# Calibration of SimpactCyan with epi-behavioural and phylogenetic tree summary statistics

Some HIV related measurements, may be elusive or difficult to estimate directly from data, thus, we have to use models, and fit those models to summary statistics (or features) from the data. In this simulation study, we investigated if some determinants of Human Immunodeficiency Virus (HIV) transmission network: (i) age-mixing patterns in sexual partnership, (ii) onward HIV transmission, and (iii) temporal trend of HIV incidence which can be computed from [SimpactCyan](<https://github.com/j0r1/RSimpactCyan/blob/master/INSTALLATION.md>), can be improved if we combine epidemiological and behavioural summary statistics with those from phylogenetic tree data.

With initial parameters’ values, we generated benchmark data, and computed values of the outcomes of interest: the determinants of HIV transmission network, and benchmark summary features from epidemiological, behavioural, and phylogenetic tree data. We considered three calibration scenarios: with benchmark features from epidemiological and behavioural data, with features from phylogenetic tree data, and with a combination of the two previous types of features. The post-calibration simulation was done by running the simulation framework using retained parameters’ values from calibration. The best calibration scenario was determined by comparing the output values after calibration to those obtained with the benchmark data through the computation of relative error.


## CONTENTS

This repo contains the information necessary to reproduce the simulation study:

* [Code and data files](#code-and-data-files)
   * Code files for simulation and post-simulation analysis
   * Data files 
   * Results files
* [System and software requirements](#system-and-software-requirements)
* [Copyright and licensing information](#copyright-and-licensing-information)
* [Contact information](#contact-information)

## CODE AND DATA FILES 

Three types of files (not including this readme file) are stored on this repository: code files, data files and figure files.


### Code files

All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

  ***benchmark*** -- contains R scripts to generate the benchmark data + one pbs file to run the simulation on CHPC
  
  ***calibration_epi_behav*** -- contains R scripts for calibration with epidemiological and behavioural summary features + one pbs file to run the simulation on CHPC
  
  ***calibration_phylo*** -- contains R scripts for calibration with phylogenetic tree summary features + one pbs file to run the simulation on CHPC

  ***calibration_combined*** -- contains R scripts for calibration with combined features from epidemiological and behavioural, and phylogenetic tree data + one pbs file to run the simulation on CHPC
  
  ***analysis_post_sim*** -- contains R scripts for data analysis and visualisation
 
 
 
 ### Data files
  
  ***data*** -- contains two subfolder: (i) ***input_data*** which contain Fasta file of input sequence for viral evolution simulation, and (ii) ***sim_outputs*** which contains sub-sub-folder for ouputs of each simulation (`benchmark`, `calibration_epi_behav`, `calibration_phylo`, and `calibration_combined`). Note that, mean vlues of summary features from benchmark data are inputs for calibration.
  
  
### Results files

  ***results*** -- tables and figures for above mentioned determinants of HIV transmission network; relative error between benchmark values and those obtained after calibration scenarios;  summary of retained parameters' values during calibration; relative error between between parameters which generated the benchmark data and those obtained with calibration; summary features associated to retained parameters' values in calibration; and relative error between benchmark summary features and those obtained with calibration.
  
  

## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system


  We run the simulation study at the University of Cape Town Centre for High Performance Computing (CHPC) and tested it on personal computer (Linux Ubuntu Version 16.04).

### Required software

  **R version 3.4.4** <www.r-project.org> For statistical computing. To Install R, do the following:
  
  1. Open an internet browser and go to www.r-project.org.
  2.  Click the "download R" link in the middle of the page under "Getting Started."
  3. Select a CRAN location (a mirror site) and click the corresponding link.
  4. Click on the "Download R for ***your OS***" link at the top of the page.
  
  

  **Seq-Gen version 1.3.4** <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network and produces a sequence alignment. To install Seq-Gen, do the following:
  
  1. Visit the following Github repository to download the latest version of Seq-Gen: <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4>
  2. Click on the "Source Code" zip file to download
  3. Click on the zip file to unzip the folder
  4. Navigate to the source folder to confirm there is a file called "Makefile"
  5. Now you will need to compile the program using the Terminal on your computer
  6. Via the Terminal, change your working directory to the source folder by typing after the prompt: `cd "file/path/here/Seq-Gen-1.3.4 2/source"`
  7. Once your working directory has been set to the source folder, type after the prompt: `make`
  8. Now open the source folder and verify that a new file is present called "seq-gen"
  9. Copy that file and paste it into your R working directory
  
  If installing on HPC facility, you may follow the instructions from 1 up to 5. And you will load the tool via the the PBS file, for example `module add /apps/chpc/scripts/modules/bio/app/Seq-Gen/1.3.4`.
  

  **FastTree version 2.1.10** <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset. To install FastTree, do the following:
  
  1. Visit the website for downloading instructions: <http://www.microbesonline.org/fasttree/#Install>
  2. If you have a Linux operating system, you can directly download the executable files that are linked on that website. Those downloaded files can then be placed in your R working directory
  3. If you are using an OS X operating system, open the link "FastTree.c" in a new browser window
  4. Right-click on the program and click "Save as"
  5. Save anywhere on your computer
  6. Open the Terminal on your computer and change your working directory to the folder that contains "FastTree.c". After the prompt type:  `cd "file/path/here"`
  7. After the directory has been changed, after the prompt type: `gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
  8. Now check to see if a new executable file has been created in that folder
  9. Copy that file and paste it into your R working directory

If installing on HPC facility, you may follow the instructions from 1 up to 7. And you will load the tool via the the PBS file, for example `module add /apps/chpc/scripts/modules/bio/app/FastTree/2.1.10`.


  **SimpactCyan version 0.21** and RSimpactCyan. SimpactCyan is the core program that allows fast simulation of HIV transmission across a sexual network. RSimpactCyan is the R package that enables initiation and running of models built by SimpactCyan. Installation instructions for both are at: <https://github.com/j0r1/RSimpactCyan/blob/master/INSTALLATION.md>


    
  A long list of auxiliary R packages is required to run the post-simulation analysis for the simulation study.

  install.packages("devtools")
  
  install.packages("pacman")
  
  library(devtools)

  install_github("j0r1/readcsvcolumns/pkg")

  install_github("wdelva/RSimpactHelp", dependencies = TRUE)

  p_load(RSimpactCyan, RSimpactHelper, Rcpp, ape, expoTree, data.table, readr, phangorn, lme4, nlme, dplyr, adephylo, treedater, geiger, picante, igraph, phyloTop, phytools, Rsamtools, robustbase, intergraph, lubridate, tidyr)
  
To run the simulation study, you need to run generate first benchmark data, by running `wrapper.benchmark.master.model.R`script, but you have to ensure that you add the executable tools "Seq-Gen", and "FastTree", to your working directory, as well as the root viral gene sequence (`hiv.seq.C.pol.j.fasta`). If you are running the simulation from HPC, you will execute `run_benchmark_master_model.pbs` file, and make sure directory is well set in the pbs file. After producing the benchmark data, we may proceed with any calibration we want maong the three scenarios. For calibration with epidemiological and behavioural summary features, you have to run `calibration_epi_behav.R` on local PC or `run_calibration_epi_behav.pbs` on CHPC. For calibration with phylogenetic tree summary features, you have to run `calibration_phylo.R` on local PC or `run_calibration_phylo.pbs` on HPC. For calibration with combined summary features, you have to run `calibration_combined.R` on local PC or `run_calibration_combined.pbs` on HPC.


**Note:** If you are running the simulation on your PC, make sure executable files of `Seq-gen` and `FastTree` are in your working directory. You need to verify also if the working directory is well set in all R files for both simulation plateforms (personal computer or HPC), and also in the pbs file if you are on HPC.




## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0 License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this work is suitable for commercial use, that licensees can modify the work, that they must release the source alongside with Derivative Work, and that Derivative Work must be released under the same terms.


## CONTACT INFORMATION

David Niyukuri
Email: <niyukuri@sun.ac.za>



