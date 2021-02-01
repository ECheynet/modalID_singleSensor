# Operational modal analysis with single accelerometer


This repository includes a Matlab function to identify the eigenfrequencies and modal damping ratios of a line-like structure using only one accelerometer. Because a single sensor is used, it is not possible to identify the mode shapes. All the data used herein are simulated.

The method applied is described as follows:

(1) The first estimate of the eigenfrequencies is obtained by the classical peak-picking method.

(2) A band-pass filter is applied to extract the modal response associated with the selected eigenfrequencies

(3) The impulse-response function (IRF) is obtained by computing the autocorrelation function of each modal response.

(4) An exponential decay function is fitted to the IRF to estimate the damping ratios and obtain an improved estimate of the eigenfrequencies, as the peak-picking method is not always appropriate.

The present approach should not be applied to structure with closely-spaced modes, as the band-pas filter may not be able to split efficiently the two modal responses.

A similar algorithm was used in [1] to identify the modal damping ratios and eigen-frequencies from traffic-induced vibrations

Reference

[1] Cheynet, E., Daniotti, N., Jakobsen, J. B., & Snæbjörnsson, J. (2020). Improved long‐span bridge modeling using data‐driven identification of vehicle‐induced vibrations. Structural Control and Health Monitoring, 27(9), e2574.
