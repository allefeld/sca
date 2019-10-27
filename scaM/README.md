# Synchronization Cluster Analysis via Markov Coarse-Graining

Code for: C. Allefeld and S. Bialonski, Detecting synchronization clusters
in multivariate time series via coarse-graining of Markov chains, *Physical
Review E*, 76:066207, 2007.

The main function is `scaM.m`. Also included are two helper functions that are
called by `scaM`, but which may also be useful for other purposes. `kMeans.m` is
a simple implementation of the k-means algorithm, `simplexGuess.m` seeks the
extreme points of a set of data points. The calling interfaces as well as
algorithmic steps are documented within the source files.

The code has been developed and tested with MATLAB and GNU Octave. If you have any problems running the code, please tell me about it.

