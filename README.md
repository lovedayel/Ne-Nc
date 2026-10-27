## Ne-Nc
__This repository contains the data and code required to recreate the figures and analysis in the manuscript: _Estimates of long-term effective population size reveal the rapid decline of endangered species___

- 'figure1.R': R code to produce main Figure 1; requires Nc_final.csv and phylo_final.nwk from the clean_data directory
- 'figure2.R': R code to produce main Figure 2; based on the analytical appendix found in the Supplementary Materials
- 'figure3.R': R code to produce main Figure 3; requires Nc_final.csv and phylo_final.nwk from the clean_data directory
- 'clean_data' directory: contains 'Nc_final.csv' (Ne, Nc & generation time estimated + IUCN status for all species included in the analysis) and 'phylo_final.nwk' (phylogenetic tree of all species included in the analysis)
- 'simulations' directory: contains the code required to reproduce the data underlying Supplementary figure 1: __run_sims.py__ (main Python script to execute simulations including the msprime burn in, reading the tree sequences into SLiM, and running the SLiM commands found in decline_instant.slim) and 'decline_instant.slim' (SLiM script called by run_sims.py to enact instantaneous population size declines and sample diversity in the subsequent generations) 
