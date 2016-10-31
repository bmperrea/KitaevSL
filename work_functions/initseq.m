function [gamma_0,gamma_pos,gamma_dec,gamma_con,var_pos,var_dec,var_con] ...
    = initseq(x);
% initseq Help file.
% [gamma_0,gamma_pos,gamma_dec,gamma_con,var_pos,var_dec,var_con] = ...
% = initseq(x);
%
% Finds the Geyer statstics of the input vector. Typically used for Monte
% Carlo algorithms.
%
% Input: vector x
% Output: [gamma_0,gamma_pos,gamma_dec,gamma_con,var_pos,var_dec,var_con]
% gamma_0 is the usual variance, var_* are better estimates of the variance 
% for a Markov chain while other gamma_* are intermediary variables in the
% calculation.
%
% See Geyer (1992) and the R function off of which this is based:
% http://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html
%
% For the actual code, see initseq.c
% To compile it, call "mex initseq.c"