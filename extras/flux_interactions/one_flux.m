rmax = 15;
s = 0.00;
b = 10;
flag = -1; %(don't plot)

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes

[H0,Rxx0,Rxy0,RxyA0,RxyB0,E0] = stretch_2D_makeH(rmax,b,s,0,0,flag);

k = 0; %The spectrum should be symmetric on k (three-fold rotation)

Es0 = zeros(rmax+1,rmax);

H = flip_plaquette(0,0,1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0); 
Es0(1,1) = -sum(svd(H))-E0;

for r = 0:rmax
    for m = 1:r
        H = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
        Es0(r+1,m) = -sum(svd(H))-E0;
    end
end

vecr = 1:(rmax+1);
vecm = 1:(rmax);
figure;
pcolor(vecr,vecm,Es0.')
colorbar
colormap hot


k = 1; %The spectrum should be symmetric on k (three-fold rotation)

Es1 = zeros(rmax+1,rmax);

H = flip_plaquette(0,0,1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0); 
Es1(1,1) = -sum(svd(H))-E0;

for r = 0:rmax
    for m = 1:r
        H = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
        Es1(r+1,m) = -sum(svd(H))-E0;
    end
end

vecr = 1:(rmax+1);
vecm = 1:(rmax);
figure;
pcolor(vecr,vecm,Es1.')
colorbar
colormap hot

Es00 = Es0;
Es11 = Es1;

