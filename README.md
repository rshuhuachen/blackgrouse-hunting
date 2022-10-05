# blackgrouse-hunting

## Workflow
Workflow for the manuscript "Effects of hunting on genetic diversity, inbreeding and dispersal in Finnish black grouse (Lyrurus tetrix)" - in prep

In this repository, you will find the scripts used for the data analysis conducted in the above manuscript.

1. Testing for Hardy-Weinberg equilibrium and removing 2 loci 
2. Analyzing population structure and calculating genetic diversity measures as well as pairwise Fst values
3. Conducting a STRUCTURE analysis and plotting the output
4. Calculating sMLH and model the effect of hunting on sMLH
5. Calculating and analyzing migration rates with BayesAss and subsequently modelling the effect of hunting on emigration and immigration rates
6. Plotting the figures from the manuscript

## Data and scripts
All data can be found within the /data directory, where raw and filtered genotypes can be found in different formats (.stru format and unsplit .csv files). In the /analysis directory you will find seperate directories for the FSTAT analysis (raw data and output), all analyses conducted in GenAlEx, an Excel add-in (used for spatial autocorrelation, sex-bias test and population assignment - not in manuscript), migration analysis (output from BA3 (BayesAss)) as well as STRUCTURE results. All R scripts can be found in the analyses/R directory, where scripts were created to be followed in order.

Parts of scripts 1 and 3 are based on the workflow for "The genetic legacy of extreme exploitation in a polar vertebrate", Paijmans et al. Sci Rep 10, 5089 (2020) which can be found at https://github.com/apaijmans/AFS_genetic_legacy.

## Checkout repository
To execute the full pipeline and copy the full repository, execute:

```
git clone https://github.com/rshuhuachen/blackgrouse-hunting.git
```

through the terminal in your own directory. Open the .Rproject file and any desired scripts. Ensure you install all required R packages as well as the GenAlEx Excel add-in, STRUCTURE and BayesAss and their dependencies.

## Additional notes

Note that in all raw genotype files, an additional locus (BG10, locus 1) is included which is not mentioned in the manuscript. We excluded this locus in R script 1 as we identified erogenous clustering in the PCA plot for females only. Most likely, a genotyping and/or scoring error occured so to avoid any biases, this locus was excluded for all analyses as well as in the manuscript.

