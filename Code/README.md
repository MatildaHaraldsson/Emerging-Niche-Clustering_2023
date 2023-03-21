# Codes

Main scripts used in simulations and calculation of clustering metrics:

* ```Main_Emerging_niche_clustering.R``` Main script to run the simulation model for the three cases: 
     - Predator-prey model with correlated competition and predator interaction niche axes (returns ```res_corr```), 
     - Predator-prey model with uncorrelated competition and predator interaction niche axes (returns ```res_uncorr```), and 
     - Competition only model without predators (returns ```res_comp```).

Using the resulting abundances for one of the simulation outputs at the time, the script runs the ```KmeansGapMod()``` algorithm to identify clusters in the trait distribution and calculates richness variables. Returns a vector of results ```SumRes``` (all successful results are found in ```datCorr.R```, ```datUncorr.R``` and ```datComp.R``` under the ```Data``` folder). It also plots the abundances of species along the niche axes with color coded clusters.

Requires the parameter values in ```parameters.R``` under the ```Data``` folder, and the functions in ```Library_Emerging_niche_clustering.R``` and ```KmeansGapMod.R```.

* ```Library_Emerging_niche_clustering.R``` Library containing functions called in the main script ```Main_Emerging_niche_clustering.R```.
* ```KmeansGapMod.R``` Contains the function ```KmeansGapMod()``` that is called in the main script  ```Main_Emerging_niche_clustering.R```. The ```KmeansGapMod()``` function is a modified version of the ```KmeansGap()``` function originally created by D’Andrea et al (2019). The code for the original metrics can be found at https://github.com/rafaeldandrea/Clustering-metric.

## References
D’Andrea R, Riolo M, Ostling A.M. (2019) Generalizing clusters of similar species as a signature of coexistence under competition. PLoS Comp Biol 15(1): e1006688
