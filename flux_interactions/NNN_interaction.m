rmax = 15;
s = 0.00;
b = 10;
flag = -1; %(don't plot)
nmax = 1;
ns = 3*nmax*(nmax+1) + 1;

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes

%Get the base Hamiltonian
[H0,Rxx0,Rxy0,RxyA0,RxyB0,~] = stretch_2D_makeH(rmax,b,s,0,0,flag);
En00 = -sum(svd(H0));

%Set a background flux.
% r=12;
% k=0;
% m=7;
% H0 = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);

%E0 = -sum(svd(H0));

numpl = rmax*(rmax+1)/2 + 1;
          %kk, l0    , l1
ps = zeros(6,numpl,ns);
rs = zeros(6,numpl,ns);
ks = zeros(6,numpl,ns);
ms = zeros(6,numpl,ns);
Ens2= zeros(6,numpl,ns);
Ens1= zeros(6,numpl);

rs0 = zeros(numpl,1);
ms0 = zeros(numpl,1);
ps0 = zeros(numpl,1);

% Single flux energies
jj = 2;
for r=1:rmax    
    for m = 1:r        
        
        rs0(jj) = r;
        ms0(jj) = m;
        ps0(jj) = 3*r*(r-1) + m + 1;
                
        for kk = 1:6 %indices for 0,1 resp.
            H1 = flip_plaquette(r,kk-1,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
            En1 = -sum(svd(H1))-En00;
            Ens1(kk,jj) = En1;
        end
        jj=jj+1;
    end
end


% center case
rs0(1)=0;
ms0(1)=0;
ps0(1)=1;
[ps1,rs1,ks1,ms1] = get_NNNNs(0,0,0,nmax);
ps(1,1,:) = ps1;
rs(1,1,:) = rs1;
ks(1,1,:) = ks1;
ms(1,1,:) = ms1;

H1 = flip_plaquette(0,0,0,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
En1= -sum(svd(H1))-En00;
for q = 1:6
    Ens1(q,1) = En1;
end       %Ens1(2,1) = En1;

for ii = 1:ns
    j0 = rs1(ii)*(rs1(ii)-1)/2 + ms1(ii) + 1;
    En0 = Ens1( mod(ks1(ii),2)+1, j0);    
    
    H = flip_plaquette(rs1(ii),ks1(ii),ms1(ii),rmax,H1,...
                                Rxx0,Rxy0,RxyA0,RxyB0);
    En = -sum(svd(H))-En00 - En1 -En0;
    Ens2(1,1,ii) = En;
end

jj=2;
%Loop over plaquettes and their NNNs
tic
for r=1:rmax
    disp(r)
    toc
    for m = 1:r
                
        for kk = 1:6 %indices for 0,1 resp.
            H1 = flip_plaquette(r,kk-1,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
            En1 = Ens1(kk,jj); %-sum(svd(H1))-En00;            
            
            [ps1,rs1,ks1,ms1] = get_NNNNs(r,kk-1,m,nmax);
            ps(kk,jj,:) = ps1;
            rs(kk,jj,:) = rs1;
            ks(kk,jj,:) = ks1;
            ms(kk,jj,:) = ms1;
            
            for ii = 1:ns
                
                if rs1(ii) <= rmax                    
                    j0 = rs1(ii)*(rs1(ii)-1)/2 + ms1(ii) + 1;
                    En0 = Ens1( mod(ks1(ii),2)+1, j0);
                    
                    %H = flip_plaquette(rs1(ii),ks1(ii),ms1(ii),rmax,H0,...
                    %                            Rxx0,Rxy0,RxyA0,RxyB0);
                    %disp(En0)
                    %disp(-sum(svd(H))-En00);
                    %En0=-sum(svd(H))-En00;
                    %Ens1(kk,jj) = En0;
                    
                    H = flip_plaquette(rs1(ii),ks1(ii),ms1(ii),rmax,H1,...
                                                Rxx0,Rxy0,RxyA0,RxyB0); 
                    En = -sum(svd(H))-En00 - En0 - En1;
                    Ens2(kk,jj,ii) = En;
                end
                
            end
        end
        jj=jj+1;
    end
end

%Save this data.
save(['flux_interaction_data_n_',num2str(nmax),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)],...
    'En00','Ens1','Ens2','ps','rs','ks','ms','rs0','ms0','ps0')

toc;

% vecr = 1:(rmax+1);
% vecm = 1:(rmax);
% figure;
% pcolor(vecr,vecm,(Es0-Es00).')
% colorbar
% colormap hot


% k = 1; %The spectrum should be symmetric on k (three-fold rotation)
% 
% Es1 = zeros(rmax+1,rmax);
% 
% H = flip_plaquette(0,0,1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0); 
% Es1(1,1) = -sum(svd(H))-E0;
% 
% for r = 0:rmax
%     for m = 1:r
%         H = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
%         Es1(r+1,m) = -sum(svd(H))-E0;
%     end
% end
% 
% vecr = 1:(rmax+1);
% vecm = 1:(rmax);
% figure;
% pcolor(vecr,vecm,(Es1-Es11).')
% colorbar
% colormap hot

