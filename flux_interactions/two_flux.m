rmax = 15;
s = 0.04;
b = 10;
flag = -1; %(don't plot)

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes

[H0,Rxx0,Rxy0,RxyA0,RxyB0,~] = stretch_2D_makeH(rmax,b,s,0,0,flag);

%Set a background flux.
r=5;
k=0;
m=0;
H0 = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
%E0 = -sum(svd(H0));

%Set a background flux pattern away from the sampling area.
% pnum = 100;
% for j = 1:pnum
%     r=randi([3,rmax]);
%     k=randi([2,4]);
%     m=randi(r);
%     H0 = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
% end
E0 = -sum(svd(H0));

%compute the spectra of fluxes around it and plot the difference from no
%flux
k = 0; %The spectrum should be symmetric on k (three-fold rotation)

Es0 = zeros(rmax+1,rmax);

H = flip_plaquette(0,0,1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0); 
Es0(1,1) = -sum(svd(H))-E0;

for r = 0:rmax
    for m = 1:r
        if r~=12 || m ~= 7
        H = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
        Es0(r+1,m) = -sum(svd(H))-E0;
        end
    end
end

Ens02 = (Es0-Es00).';

save(['flux_interaction_data1_inf_',num2str(nmax),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)],...
    'E0','Es0','Ens02','Es00','r','k','m')

vecr = 1:(rmax+1);
vecm = 1:(rmax);
figure;
pcolor(vecr,vecm,(Es0-Es00).')
colorbar
colormap hot

hh=figure;
hexagonal_density_plot(rmax,s,Ehist)
filename = ['E2hist_Tc_',num2str(round(10000*Tv(minT))),'_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename); 


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

