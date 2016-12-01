function Z = getPartitionFunction( dos, Ev, Tv )
%getPartitionFunction Gets the partition function as a function of
    %temperature as an afterthought to flux_controller_final.
%   sites is the number of sites
%   dos   is the expected DOS as a function of T,E
%   Ev    is the energy values for the corresponding energy asix in DOSs
%   Tv    is the vector of temps

%Typical usage: ZZ = getPartitionFunction( II{2}, Ev, Tv );

% Based on: Z = sum_p BF_p sum_w (1+exp(-beta w)) <dos_p>(w)
           %  = sum_w (1+exp(-beta w)) sum_p <dos_p>(w)
           %  = sum_w (1+exp(-beta w)) <dos>(w).

    Z = zeros(numel(Tv),1);
    for in = 1: numel(Z)
        bf0 = reshape( (1+exp(-(Ev/2) / Tv(in))) , 1,numel(Ev) );
        dE = max(Ev/2)/200;
        Z(in) = dE * sum( bf0 .* dos(in,:) ) ;
        Z(in) = Z(in);
    end

end