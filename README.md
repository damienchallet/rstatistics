# rstatistics
Code for computing r-statistics, a new nonparametric statistics based on the number of records of the cumulative sum of sample points. 

The code allows you to use straightforwardly this new statistics. For the time being, i.e., until I create an R package, to compute r-statistics, type

library(Rcpp)

sourceCpp('num_records.cpp')

For single samples, rstat(x) returns the r-statistics of vector x. 

For two-sample situations (say x and y), there are several statistics. The simplest idea is to use rstat(x-y). Otherwise, computeUpsDownsTwoSamplesRandomPerms(x,y) returns several quantities from which one can compute r-statistics (TODO: make it more user-friendly).

This code also makes it possible reproduce results of two submitted papers of mine, http://arxiv.org/abs/1505.01333 and http://arxiv.org/abs/1502.05367, including their plots.
