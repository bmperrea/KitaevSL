function En = predict_En0(r,kk,m,~,Ens1,~,~)

%    kk = mod(k,2);
    j = r.*(r-1)/2 + m + 1;
    
    En1 = zeros(size(j));
    
    %if pbool(j)   
    for a = 1:numel(j)
        En1(a) = Ens1(kk(a)+1,j(a));
    end
    %En1(~pbool(j)) = -En1(~pbool(j));
    %else
        %En1 = Ens1(kk+1,j);
    %end
    %Note that if there is already a flux here it will be taken care of by
    %the flux-flux interaction En2
    
%    ps = get_NNNs(r,k,m);
%    ps = ps(ps<=numplaq); %remove the stuff outside the lattice

%    En2 = sum( Ens2(kk+1,j,pbool(ps)) );
    
    En = En1;% + En2;    

end