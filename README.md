# Scripts for Novel metrics of spatial co-occurrence paper
## Description
Scripts belonging to the paper titled Novel metrics of spatial co-occurrence reflect species interaction types in a cold-water coral reef community, by A. van der Kaaden et al. 

The script Scores_Simulations performs 100 repetitions of the Lotka-Volterra style simulations from section 2.3 of the paper. It calculates the Aggregation and Symmetry Scores, calculates the quantiles of the results, and plots the results (figure 2).

The script ScoresForAnnotatedImages performs the spatial analysis (section 3.3.1) and co-occurrence plot analysis (section 3.3.2) of the images annotated by Maier et al. (2021 - see paper for full reference). It plots figure 4 and 5.

The script Scores_Matrices draws matrices with idealized distributions of two species and calculates Aggregation and Symmetry Scores for different resolutions (section 2.2 and figure 1). The script also calculates the C-score (supplement 1). 

## Installation

The code files can be downloaded and executed using R. The scripts were written and executed in R version 4.4.2 and R studio version 2025.09.2+418.  

### Dependencies

The scripts have the following dependencies:
ReacTran (Score_Simulations.R)
RColorBrewer (Score_Simulations.R, ScoresForAnnotatedImages.R)
ggplot2 (all)
png (ScoresForAnnotatedImages.R)
readxl (ScoresForAnnotatedImages.R)
terra (all)

These R packages can be installed through cran with the command install.packages(), e.g., install.packages('ReacTran').

## Usage

The scripts Score_Simulations.R and Scores_Matrices.R can be executed as they are. 
Score_Simulations.R contain a function (CalXrab) that calculates the components of the metrics across spatial scales. This function can be copied and used in new R scripts to easily calculate the metrics. 
The script ScoresForAnnotatedImages.R requires that the working directory is set to a folder with images. For the purpose of using the scripts, seven annotated images are published in the folder "images".
These are the same images as in supplementary figure 7 of the paper. 

## Citation

For the correct citations, we refer the reader to the Zenodo entry.

## How to contribute

Feel free to contact us for questions, advise, collaborations, or feedback. 
dr. Anna van der Kaaden, a.vanderkaaden@uu.nl, OrcID: 0000-0002-8814-1822
