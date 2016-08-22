# FHR_classification

============= Calculating Features =============

For calculating features, run the script 'testCompute.m'.

Note that a data file 'fhr.mat' must exist. This file should contain an array named 'data'. 
The format of 'data' should be [M x (time length)].

This script will calculate 2 minute long observations, with each observation (after the first) containing
one minute of overlap with the previous observation.

Note that I did not program error checking into this, so if a parameter is set to an invalid value it might
not throw an error but erronous data may result.

============= Training Models =============

To train a HMM model, run the script 'crossValHMM.m'.

NOTE: Training an HMM depends on the library written my Kevin Murphy that can be downloaded here:
https://www.cs.ubc.ca/~murphyk/Software/HMM/hmm.html

To train a Naive Bayes model, run the script 'crossValNB.m'

M = number of fetuses
D = number of features
N = number of observations

NOTE: The Naive Bayes script depends on a function that was not written by me. It is the multivariate log-normal pdf.
For obvious purposes I am not including this file in my upload, but it can be downloaded from here:
http://www.mathworks.com/matlabcentral/fileexchange/34064-log-multivariate-normal-distribution-function

