# Synchronization Cluster Analysis

Code for: C. Allefeld and J. Kurths, An approach to multivariate phase
synchronization analysis and its application to event-related potentials,
*International Journal of Bifurcation and Chaos*, 14(2):417-426, 2004.

The code of the SCA algorithm is `sca.m`. It takes as input the matrix of
synchronization strengths, which you have to calculate from your data before.
The interface of the function is commented in the file, and there are more
comments in the code.

This original SCA algorithm assumes that there is only one cluster in the data.
A newer approach tries to generalize to the case of several coexisting
clusters. If this is interesting to you, you should have a look at the paper
"Detecting synchronization clusters in multivariate time series via
coarse-graining of Markov chains".

The code has been developed and tested with MATLAB and GNU Octave. If you have any problems running the code, please tell me about it.

