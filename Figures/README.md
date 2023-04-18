# Figures

## Figure Code
Codes to reproduce figures 1-4 in the main text of Haraldsson and Thébault 2023. The codes demands the data objects below, as well as the ```Library_Emerging_niche_clustering.R``` script stored in the ```Code/``` folder, and the output stored in the ```Data/``` folder.

- ```Figure_1.R```: reproduces figure 1
- ```Figure_2_3.R```: reproduces figure 2 and 3
- ```Figure_4.R```: reproduces figure 4

## Figure Data
- ```datClust_spnr2_t500_852226.RData```: clustering results from the simulation (simulation ID 852226) as show in Figure 1.
- ```nbclustmat_Q4_aR0_01.RData```: expected number of clusters as function of sigma N and sigma P, with a "Importance of Predation" = 0.01, as shown in Figure 1.
- ```nbclustmat_Q4_aR1.RData```: expected number of clusters as function of sigma N and sigma P at a "Importance of Predation" = 1, as shown in Figure 1.
- ```nbclustmat_Q4_aR1.RData```: expected number of clusters as function of sigma N and sigma P at a "Importance of Predation" = 100, as shown in Figure 1.
- ```long_abu_comp.RData```: richness of species with time in the competing community, as shown in Figure 4.
- ```long_abu_corr_N.RData```: richness of species with time in the competing prey community in the correlated case, as shown in Figure 4.
- ```long_abu_corr_P.RData```: richness of species with time in the predator community in the correlated case, as shown in Figure 4.
- ```long_abu_uncorr_N.RData```: richness of species with time in the competing prey community in the uncorrelated case, as shown in Figure 4.
- ```long_abu_uncorr_NP.RData```: richness of species with time in the predator community in the uncorrelated case, as shown in Figure 4.
- ```long_clust_comp_shuffle.RData```: clustering results with time in the competing community, using the "shuffle" method, as shown in Figure 4.
- ```long_clust_corr_shuffle.RData```: clustering results with time in the competing prey community in the correlated case, using the "shuffle" method, as shown in Figure 4.
- ```long_clust_uncorr_shuffle.RData```: clustering results with time in the predator community in the correlated case, using the "shuffle" method, as shown in Figure 4.
