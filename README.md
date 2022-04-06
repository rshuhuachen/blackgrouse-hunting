# blackgrouse-hunting
Workflow for "Sex-specific fine-scale population structure and effects of hunting on the genetic diversity and dispersal of Finnish black grouse (Lyrurus tetrix)" - in prep

In this repository, you will find the 5 scripts used for the data analysis conducted in the above manuscript.

1. Testing for Hardy-Weinberg equilibrium in adult data
2. Analyzing population structure through making a PCA and calculating pairwise Fst values
3. Conducting a STRUCTURE analysis and plotting the output
4. Calculating sMLH for adults and chicks and model the effect of hunting on sMLH
5. Analyzing migration rates by modelling the effect of hunting on emigration and immigration
6. Plotting the figures from the manuscript

All data can be found within the /data directory, where raw and filtered genotypes can be found in various formats (stru format and unsplit csv files), pairwise Fst values, spatial autocorrelation autput, migration analysis output from BA3 (Bayesass), STRUCTURE results (partly), and sMLH per individual. To get a better description for each data file, please follow the scripts 1-5 in order.

Parts of scripts 1 and 3 are based on the workflow for "The genetic legacy of extreme exploitation in a polar vertebrate", Paijmans et al. Sci Rep 10, 5089 (2020) which can be found at https://github.com/apaijmans/AFS_genetic_legacy.
