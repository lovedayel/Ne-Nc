## Ne-Nc
# This repository contains the data and code required to recreate the figures and analysis in the manuscript: _Estimates of long-term effective population size reveal the rapid decline of endangered species_ 

- __figure1.R__: R code to produce main Figure 1; requires Nc_final.csv and phylo_final.nwk from the clean_data directory
- __figure2.R__: R code to produce main Figure 2; based on the analytical appendix found in the Supplementary Materials
- __figure3.R__: R code to produce main Figure 3; requires Nc_final.csv and phylo_final.nwk from the clean_data directory
- __clean_data__ directory: contains __Nc_final.csv__ (Ne, Nc & generation time estimated + IUCN status for all species included in the analysis) and __phylo_final.nwk__ (phylogenetic tree of all species included in the analysis)
- __simulations__ directory: contains the code required to reproduce the data underlying Supplementary figure 1: __run_sims.py__ (main Python script to execute simulations including the msprime burn in, reading the tree sequences into SLiM, and running the SLiM commands found in decline_instant.slim) and __decline_instant.slim__ (SLiM script called by run_sims.py to enact instantaneous population size declines and sample diversity in the subsequent generations) 
