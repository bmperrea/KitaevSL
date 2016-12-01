
hh=figure;
pp = (1:850)*.5/1000;
tcs = -1./log2( pp./(1-pp) );
plot(pp,tcs)
xlabel('p_c');
ylabel('T_c/\Delta');
filename = ['Tc_toy_model'];
axis([0,.5,0,inf])
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);  


hh=figure;
%axis square
pp = (1:1000)/1000;
tcs = -pp.*log2( pp) - (1-pp).*log2(1-pp);
plot(pp,tcs)
xlabel('\rho');
ylabel('S/N_P');
filename = ['S_toy_model'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);  

s = 0.00; rmax = 15; numplaq = 1 + 3*rmax*(rmax+1); b =10;

[rus,kus,mus] = get_polar(reshape(1:numplaq,numplaq,1));
if s == 0.04
    nmax = 5;
else
    nmax =1;
end
load(['flux_interaction_data_n_',num2str(nmax),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)])
Ehist = predict_En0(rus,mod(kus,2),mus,zeros(numplaq,1),Ens1);

hh=figure;
hexagonal_density_plot(rmax,s,Ehist)
filename = ['E1hist_Tc_',num2str(round(10000*Tv(minT))),'_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);   

% Ehist2 = zeros(1,numplaq);
% for a = 1:numplaq
%     [ru,ku,mu]=get_polar(a);
%     Ehist2 = predict_En(ru,mod(ku,2),mu,pbool,Ens1,Ens2,5);
% end
% 
% hh=figure;
% hexagonal_density_plot(rmax,s,Ehist)
% filename = ['E1hist_Tc_',num2str(round(10000*Tv(minT))),'_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);   




