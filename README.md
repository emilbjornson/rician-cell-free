Performance of Cell-Free Massive MIMO with Rician Fading and Phase Shifts
==================

This is a code package is related to the follow scientific article:

Özgecan Özdogan, Emil Björnson, Jiayi Zhang, “[Performance of Cell-Free Massive MIMO with Rician Fading and Phase Shifts](https://arxiv.org/pdf/1903.07335.pdf),” IEEE Transactions on Wireless Communications, To appear.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

In this paper, we study the uplink (UL) and downlink (DL) spectral efficiency (SE) of a cell-free massive multiple-input-multiple-output (MIMO) system over Rician fading channels. The phase of the line-of-sight (LoS) path is modeled as a uniformly distributed random variable to take the phase-shifts due to mobility and phase noise into account. Considering the availability of prior information at the access points (APs), the phase-aware minimum mean square error (MMSE), non-aware linear MMSE (LMMSE), and least-square (LS) estimators are derived. The MMSE estimator requires perfectly estimated phase knowledge whereas the LMMSE and LS are derived without it.

In the UL, a two-layer decoding method is investigated in order to mitigate both coherent and non-coherent interference. Closed-form UL SE expressions with phase-aware MMSE, LMMSE, and LS estimators are derived for maximum-ratio (MR) combining in the first layer and optimal large-scale fading decoding (LSFD) in the second layer. In the DL, two different transmission modes are studied: coherent and non-coherent. Closed-form DL SE expressions for both transmission modes with MR precoding are derived for the three estimators. Numerical results show that the LSFD improves the UL SE performance and coherent transmission mode performs much better than non-coherent transmission in the DL. Besides, the performance loss due to the lack of phase information depends on the pilot length and it is small when the pilot contamination is low.


## Content of Code Package

See each file for further documentation.


## Acknowledgements

This work was supported in part by ELLIIT and the Swedish Research Council. It was also supported by National Natural Science Foundation of China (Grant Nos. 61601020 and U1834210), the Beijing Natural Science Foundation (Grant Nos. 4182049 and L171005).


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
