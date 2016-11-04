function [var_con,var_dec,var_pos,G_con,G_dec,G_pos,G_0] ...
    = initseq_matlab(x)

% Finds the Geyer statstics of the input vector.
%
% Input: vector x
% Output: [gamma_0,gamma_pos,gamma_dec,gamma_con,var_pos,var_dec,var_con]
% var_* are estimates of the variance for a Markov chain 
% gamma_* are intermediary variables in the calculation.
%
% See Geyer (1992) and the R function off of which this is based:
% http://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html
%
% This version is purely written with standard Matlab functions and 
% competes well with the mex file initseq_c for performance.
%
% There is a "vectorized" version that calls this in a loop called
% initseq_matlab_vec.m 
% 
% There is a faster version that computes batch means to reduce the data
% before the computation in initseq_batch.m

    %Get sizes
    len = numel(x);
    
    %subtract means
    x = bsxfun(@minus, x, mean(x));
    
    %Find autocorrs (gamma)
    [acf,lags] = xcorr(x);
    ac = acf(lags>=0);
    G_0 = ac(1) / len;
    
    %Find Gamma
    ii = floor(len/2);
    Gamma = sum(reshape(ac(1:(2*ii)),2,ii),1) / len;
    
    %Find ipse
    inds = find(Gamma<0);
    if numel(inds) == 0
        inds = numel(Gamma)+1;
    end
    G_pos = [Gamma(1 : max(inds(1)-1, 1) ), 0];
    var_pos = 2*sum(G_pos) - G_0;    
    
    %Find imse
    G_dec = cummin(G_pos);
    var_dec = 2*sum(G_dec) - G_0;
        
    %Find icse
    buff = diff(G_dec);
    diff_pos = pav(buff);
    G_con = G_dec(1) + [0, cumsum(diff_pos)];
    var_con = 2*sum(G_con) - G_0;

end

function mincon = pav(x)
% Pool Adjacent Violators (PAV)
    mincon = lsqisotonic((1:numel(x)),x);
end