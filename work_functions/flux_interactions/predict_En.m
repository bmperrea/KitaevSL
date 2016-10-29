function [En,En0] = predict_En(r,k,m,pbool,Ens1,Ens2,numplaq,nmax)

 %   kk = mod(k,2);
    j = r*(r-1)/2 + m + 1;
    
    En1 = Ens1(k+1,j);
    %Note that if there is already a flux here it will be taken care of by
    %the flux-flux interaction En2
        
    if nargout == 2
        En0 = predict_En0(r,k,m,pbool,Ens1);
    end
    
    ps = get_NNNNs(r,k,m,nmax);
    ps = ps(ps<=numplaq); %remove the stuff outside the lattice
    flipped = pbool(ps(1));
    pbool(ps(1)) = 0;

    En2 = sum( Ens2(k+1,j,pbool(ps)) );
    
    En = En1 + En2;    
        
    if flipped
        En = -En;
    end
    
end