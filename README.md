# Operational modal analysis with a single accelerometer

[![View Operational modal analysis with a single accelerometer on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/86718-operational-modal-analysis-with-a-single-accelerometer)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4487060.svg)](https://doi.org/10.5281/zenodo.4487060)
<a href="https://www.buymeacoffee.com/echeynet" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 25px !important;width: 120px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

## Summary
This repository includes a Matlab function to identify the eigenfrequencies and modal damping ratios of a line-like structure using only one accelerometer. Because a single sensor is used, it is not possible to identify the mode shapes. All the data used herein are simulated. I am using the peak picking function “pickpeaks” developed by [1] and available in [2]. The manual pick-peaking algorithm was proposed by [3] in the Matlab file exchange.  I am, therefore, indebted to [1-3] for their previous works. 

The method applied is described as follows:

(1) The first estimate of the eigenfrequencies is obtained by the classical peak-picking method.

(2) A band-pass filter is applied to extract the modal response associated with the selected eigenfrequencies

(3) The impulse-response function (IRF) is obtained by computing the autocorrelation function of each modal response.

(4) An exponential decay function is fitted to the IRF to estimate the damping ratios and obtain an improved estimate of the eigenfrequencies, as the peak-picking method is not always appropriate.

The present approach should not be applied to structure with closely-spaced modes, as the band-pass filter may not be able to split efficiently the two modal responses.

A similar algorithm was used in [4] to identify the modal damping ratios and eigenfrequencies from traffic-induced vibrations. In [4], no advanced pick-peaking algorithm was used as we already knew what were the eigenfrequencies of the structure.


## Content

The repository contains:
  - The function modalID_sngleSensor.m, which identify the eigenfrequencies and damping ratios of a line-like structure
  - A Matlab livescript Documentation.mlx
  - Some simulated data BridgeAccData.mat for the examples in Documentation.mlx.


## References

[1] Antoine Liutkus. Scale-Space Peak Picking. [Research Report] Inria Nancy - Grand Est (Villers-lès-Nancy, France). 2015. .

[2] https://se.mathworks.com/matlabcentral/fileexchange/42927-pickpeaks-v-select-display-

[3] https://se.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-

[4] Cheynet, E., Daniotti, N., Jakobsen, J. B., & Snæbjörnsson, J. (2020). Improved long‐span bridge modeling using data‐driven identification of vehicle‐induced vibrations. Structural Control and Health Monitoring, 27(9), e2574.
