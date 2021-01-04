# Spatio-temporal food web dynamics

Romain Frelat, 4th January 2020



The repository contains the scripts and dataset to reproduce the results presented in the manuscript "Spatio-temporal community dynamics: the missing link between composition and structure" by Romain Frelat, Susanne Kortsch, Ingrid Kröncke, Hermann Neumann, Marie C. Nordström, Pierre Olivier, and Anne F. Sell.



The scripts are commented as much as possible, but are not tutorials. For tutorials about food web analysis, please visit: https://rfrelat.github.io/BalticFoodWeb.html. For tutorials about tensor decomposition, please visit: https://rfrelat.github.io/Multivariate2D3D.html.



The analysis unfold in three steps:

* 1_TensorComposition.R: Tensor decomposition on the relative abundance and clustering of the species based on their spatio-temporal dynamics
* 2_TensorStructure.R: Calculation of the food web metrics per box and per year followed by the tensor decomposition of the food web metrics.
* 3_CoinertiaPTA.R: Co-inertia of the two tensor decomposition, linking the spatio-temporal dynamics of community composition and structure.



Packages to run these scripts are: ade4, PTAk, igraph, RColorBrewer and corrplot. If some of these packages are not installed, you can install them with the following command line:

```{r}
install.packages(c("ade4", "PTAk", "igraph", "RColorBrewer", "corrplot"))
```

Version: all analyses were run on R4.0.2.



The dataset and additional functions are all included in the Rdata file [TensorNorthSea.Rdata](https://github.com/rfrelat/NorthSeaFoodWeb/raw/main/TensorNorthSea.Rdata) which contains:

- tensorNS: the tensor of relative abundance of 114 species, over 6 boxes and 15 year
- netNS: food web of the north Sea of the 114 persistant species and five functional groups (Detritus,  Microalgae, Macroalgae, Phytoplankton, Zooplankton)
- coordinatesBox: coordinates of the 6 boxes (longitude and latitude in decimal degrees)
- six home-made functions (use them at your own risk):
  - trophiclevels: calculate the trophic level
  - jaccardsim: calculate the diet similarity, defined as the Jaccard distance
  - fwind: extract multiple food web metrics
  - myHeatmap: plotting function to custumize a heatmap
  - maploadings: plotting function to plot the loadings based on coordinates





**Be aware**: the relative abundance were derived from pre-processing of the raw data to fit the purpose of our study (Hellinger log transformed and constant sampling intensity, see Material and Methods for more details). To avoid any misuse of the dataset, the species names were replaced by numbers (from S1 to S114) in tensorNS and netNS. If you require the raw dataset from the German Small Scale Bottom Trawl Survey, please contact Anne Sell and Ingrid Kröncke.





*Add DOI*

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

