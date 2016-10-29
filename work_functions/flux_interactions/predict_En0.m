function En = predict_En0(r,kk,m,pbool,Ens1,~,~)

%    kk = mod(k,2);
    j = r*(r-1)/2 + m + 1;
    
    if pbool(j)
        En1 = -Ens1(kk+1,j);
    else
        En1 = Ens1(kk+1,j);
    end
    %Note that if there is already a flux here it will be taken care of by
    %the flux-flux interaction En2
    
%    ps = get_NNNs(r,k,m);
%    ps = ps(ps<=numplaq); %remove the stuff outside the lattice

%    En2 = sum( Ens2(kk+1,j,pbool(ps)) );
    
    En = En1;% + En2;    

end