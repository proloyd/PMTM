# MultitaperSpectrum_from_BinaryNeuralData

Description: This repository contains implementations of the multitaper PSD estimation method developed for binary spiking data.

Copyright (c) 2017 Proloy Das All Rights Reserved 

Contact: proloy@umd.edu

(Das and Babadi, Multitaper Spectrum from Binary Neural Data,(under review) Submitted @ IEEE Signal Process. Lett. on )

Date: June 5, 2017

Requirements:
  implemented in Matlab R2013b version, but should run on most versions.
  
Contents:

    1. PMTM.m: the main script.
    3. SpikeRaterPlot.m: creates raster plot from Binary data.
    4. SS_latent_estimation.m: estimate latent process under AR(1) assumption, (used to generate SS-PSD).
    
Instructions: Simple and easy.
  Download all the codes in a directory and run the multitaperPsdBinarySpikes.m file, that will generate AR process, binary spiking data from it and finally the multitaper PSD estimate from the binary data. At the end, it will generate an overlay plot(similar to Fig. 3A in the paper) to compare it to other estimates.
