# PMTM Point Process Multitaper Method

Description: This repository contains implementations of the multitaper PSD estimation method developed for binary spiking data.

[1] Das, P., & Babadi, B. (2019). [Multitaper Spectral Analysis of Neuronal Spiking Activity Driven by Latent Stationary Processes.](https://arxiv.org/pdf/1906.08451) arXiv preprint arXiv:1906.08451.

Copyright (c) 2017 Proloy Das All Rights Reserved 

Contact: proloy@umd.edu

Requirements:
  implemented in Matlab R2013b version, but should run on most versions.
  
Contents:

    1. PMTM.m: the main script.
    3. SpikeRaterPlot.m: creates raster plot from Binary data.
    4. SS_latent_estimation.m: estimate latent process under AR(1) assumption, (used to generate SS-PSD).
    
Instructions: Simple and easy.
  Download all the codes in a directory and run the PMTM.m file, that will generate a set of binary spiking data and compute the multitaper PSD estimate from it. At the end, it will generate an overlay plot (Fig. 3A in the paper) to compare it to other estimates.
