rmax = 15;
s = 0.00;
b = 10;
flag = -1; %(don't plot)

if s == 0.04
    nmax = 5;
else
    nmax =1;
end
load(['flux_interaction_data_n_',num2str(nmax),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)])

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes

[H0,Rxx0,Rxy0,RxyA0,RxyB0,~] = stretch_2D_makeH(rmax,b,s,0,0,flag);
E00 = -sum(svd(H0));

%Set a background flux.
r1=5;
k1=0;
m1=1;
H0 = flip_plaquette(r1,k1,m1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
%E0 = -sum(svd(H0));

%Set a background flux pattern away from the sampling area.
% pnum = 100;
% for j = 1:pnum
%     r=randi([3,rmax]);
%     k=randi([2,4]);
%     m=randi(r);
%     H0 = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
% end
E0 = -sum(svd(H0))-E00;

%compute the spectra of fluxes around it and plot the difference from no
%flux
%k = 0; %The spectrum should be symmetric on k (three-fold rotation)

%Es0 = zeros(rmax+1,rmax);

%H = flip_plaquette(0,0,1,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0); 
%Es0(1,1) = -sum(svd(H))-E0;

Ens02 = zeros(numplaq,1);
lastr = 0;
for a = 1:numplaq
    [r,k,m]=get_polar(a);
    if lastr ~= r %This is just to see progress...
        lastr = r;
        disp(r)
    end
   % if r~=r1 || m ~= m1
        En1 = predict_En0(r,k,m,0,Ens1);
        H = flip_plaquette(r,k,m,rmax,H0,Rxx0,Rxy0,RxyA0,RxyB0);
        Ens02(a) = (-sum(svd(H))-E00) - E0 - En1;
 %   end
end

indy=3*r1*(r1-1) +1  + k1*rmax + m1;
Ens02(indy) = 0;

%Ens02 = (Es0-Es00).';

save(['flux_interaction_data1_',num2str(r1),'_',num2str(k1),'_',num2str(m1), ...
    '_inf_',num2str(nmax),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)],...
    'E0','Ens02','E00','r1','k1','m1')

% vecr = 1:(rmax+1);
% vecm = 1:(rmax);
% figure;
% pcolor(vecr,vecm,(Es0-Es00).')
% colorbar
% colormap hot

%sorty = sort(Ens02);
%minny = sorty(2);

hh=figure;
hexagonal_density_plot(rmax,s,Ens02)

sf = 3*(1:64)/64; cmap2 = min(max(0, [sf/2 ; sf*2-4 ; sf*2-5 ]),1).'; 
cmin = min(Ens02);
cmax = max(Ens02);
caxis([cmin,cmax])
mind = round( 64*cmax/(cmax-cmin) );
%sf2 = sf(1:mind)*2;%(1:mind)/mind;
%negs = find( cmap3(:,1)<0 );
if mind > 21  
    sf2 = (1:mind)/mind;
else
    sf2 = sf(1:mind)*2;%(1:mind)/mind;
end
cmap3 = [ cmap2(mind:end,:) ; 1-[sf2;sf2;sf2*0].' ];
colormap(cmap3)

%The actual positions after strain
r0 = rr(r1,k1,m1,0);     %even sublattice
r5 = rr(r1-1,k1,m1-1,1);   %odd
rav = (r0 + r5)/2;
text(rav(1)-.7,rav(2)+.4,'\pi')

filename = ['E2hist_',num2str(r1),'_',num2str(k1),'_',num2str(m1),...
    '_Tc_',num2str(round(10000*Tv(minT))),'_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
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

