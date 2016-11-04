function [var_pos,var_dec,var_con,gamma_0] = initseq_c_vec(x)
% initseq Help file.
% [var_pos,var_dec,var_con,gamma_0] = initseq_vec(x)
%
% Finds the Geyer statstics of the input vector. Typically used for Monte
% Carlo algorithms.
%
% Input: x   matrix whose columns are Markov Chains.
% Output: [var_pos,var_dec,var_con,gamma_0]
% row vectors representing variance estimates for each column of x
% var_* are estimates of the variance for a Markov chain. 
%
% See Geyer (1992) and the R function off of which this is based:
% http://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html
%
% For the actual code, see initseq_vec.c
% To compile it, call "mex initseq_vec.c"