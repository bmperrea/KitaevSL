tic

prep_plots
close all
savePlots = true;

favColor1 = [153,216,201]/255 - .25;
favColor2 = [44,152,105]/255 - .15;
lw = 2;

rmax = 15; nt = 2.5;
s = .00;% * (rmax + 1/2)/(15 + 1/2);
b=10;
numplaq = 3 * rmax * (rmax +1 ) +1;

fileName = ['stretch_flux_mcmc1_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s))];
load(fileName)    

%Tc = 2.5*Tc;

format compact
%Get flux sector distribution and energy distributions

flux_dist = zeros(numel(Tv),numplaq+1);
energy_dist = zeros(numel(Tv),numplaq+1);

if rmax > 9 
    enpplaq = 0.06;
else 
    enpplaq = 0.07;
end

maxE = 0;
maxp = 0;

pvalys = zeros(1,50000*50);
Evalys = zeros(1,50000*50);
ii = 1;

counts0 = zeros(1,numel(Tv));
errors0 = counts0;

phists = zeros(numel(Tv),numplaq);
[rus,kus,mus] = get_polar(reshape(1:numplaq,numplaq,1));
rhists = zeros(numel(Tv),rmax+1);
khists = zeros(numel(Tv),6);
mhists = zeros(numel(Tv),rmax+1);

[tmin,minT]=min(abs(Tv-Tc));
disp(Tv(minT))

for ind = 1:numel(Tv)
    
    T = Tv(ind);
    load(['stretch_mcmc_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),...
        '_s_',num2str(floor(1000*s)),'_T_',num2str(round(T*10000))])
    
    flux_dist(ind,:) = accumarray(pss(start:end)+1,1,[numplaq+1,1]).' / (count-start);
        
    enis = round(Ens(start:end) / enpplaq );
    energy_dist(ind,:) = accumarray(enis+1,1,[numplaq+1,1]).' / (count-start);
    
    maxE = max(maxE,max(enis));
    maxp = max(maxp,max(enis));
    
    pmin = min(pss(start:end)+1);
    pmax = max(pss(start:end)+1);
    
    for jj = 1:(count-start+1)
        pvalys(ii) = pss(jj);
        Evalys(ii) = Ens(jj);
        ii = ii + 1;
    end
    
    if error(4) > 0.02
        disp([ind,count])
        disp(error)
    end
    
    errors0(ind) = max(error);
    
    counts0(ind) = count;
    
    phist = phist/count;    
    rhist = accumarray(rus+1,phist)./( [1,6*(1:rmax)].' );
    khist = accumarray(kus+1,phist);
    mhist = accumarray(mus+1,phist)./( [1,6*flip(1:rmax)].' );
    
    phists(ind,:) = phist;
    rhists(ind,:) = rhist;
    khists(ind,:) = khist;
    mhists(ind,:) = mhist;
     
    if ind == minT       
        hh=figure;
        hexagonal_density_plot(rmax,s,phist)    
        %xlabel('T/J');
        %ylabel('p');
        filename = ['phist_Tc_',num2str(round(10000*Tv(minT))),'_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);    
    end
    
end

disp(counts0)
disp(errors0)

maxp2 = round(1/2 * numplaq);
maxE2 = Tc;

hh=figure;%('Position',position);
hold on;
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), (0:numplaq)/numplaq, flux_dist(Tvis,:).' );
%axis([0, inf, -inf,inf])
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, maxp2/numplaq])
m3 = 0.15;%max(max( flux_dist(Tvis,:) ));
caxis([0,m3])
xlabel('T/J');
ylabel('\rho');
hold off;

    filename = ['p_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);

hh=figure;%('Position',position);
hold on;
ensi = (0:numplaq) * enpplaq / (numplaq+1);
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), ensi, energy_dist(Tvis,:).' );
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, maxE2])
m3 = 0.15;%max(max( energy_dist(Tvis,:) ));
caxis([0,m3])
xlabel('T/J');
ylabel('E_0/J');
hold off;

    filename = ['E_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);

    %phist
hh=figure;%('Position',position);
hold on;
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), 1:numplaq, phists(Tvis,:).' );
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, numplaq])
%m3 = 0.15;%max(max( energy_dist(Tvis,:) ));
%caxis([0,m3])
xlabel('T/J');
ylabel('pnum');
hold off;

    filename = ['pnum_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);    
    
    %rhist
hh=figure;%('Position',position);
hold on;
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), 1:(rmax+1), rhists(Tvis,:).' );
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, (rmax+1)])
%m3 = 0.15;%max(max( energy_dist(Tvis,:) ));
%caxis([0,m3])
xlabel('T/J');
ylabel('r');
hold off;

    filename = ['r_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename); 
    
    %khist
hh=figure;%('Position',position);
hold on;
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), 1:6, khists(Tvis,:).' );
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, 6])
%m3 = 0.15;%max(max( energy_dist(Tvis,:) ));
%caxis([0,m3])
xlabel('T/J');
ylabel('k');
hold off;

    filename = ['k_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename); 

 %mhist
hh=figure;%('Position',position);
hold on;
Tvis = (Tv <= (nt+1)*Tc);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor(Tv(Tvis), 1:(rmax+1), mhists(Tvis,:).' );
colorbar
colormap hot
shading flat
axis([0, nt*Tc, 0, (rmax+1)])
%m3 = 0.15;%max(max( energy_dist(Tvis,:) ));
%caxis([0,m3])
xlabel('T/J');
ylabel('m');
hold off;

    filename = ['m_v_T_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);    
    
     %Energy versus p
     
%     hh=figure;
%     hold on;
%     scatter(pvalys(1:ii)/numplaq,Evalys(1:ii)/numplaq,1,'filled')
%     xlabel('\rho');
%     ylabel('Energy per plaq.');
%     hold off;
%     filename = ['stretch_En_mcmc_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(floor(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end

maxp3 = round(3/5 * numplaq);
    
    Evals = Evalys(1:ii);
    pvals = pvalys(1:ii);
    
    energy_pdist = zeros(maxp3+1,numplaq+1);
    for ind = 0:maxp3
        enis = round( Evals(pvals == ind) / enpplaq ).';
        if numel(enis) > 0
            energy_pdist(ind+1,:) = accumarray(enis+1,1,[numplaq+1,1]).' / (count-start);
        end
    end
hh=figure;%('Position',position);
hold on;
ensi = (0:numplaq) * enpplaq / (numplaq+1);
% pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
pcolor((0:maxp3)/numplaq, ensi, energy_pdist.' );
colorbar
colormap hot
shading flat
axis([0, maxp3/numplaq, 0, Tc*1.25])
m3 = max(max( energy_pdist ));
caxis([0,.01])
xlabel('\rho');
ylabel('(E_0/J)/N_p');
hold off;
filename = ['E_v_p_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);    
    
    
%rmax = 15;
%s = .04;
fileName = ['stretch_flux_mcmc1_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s))];
load(fileName)
disp(max(Ev))
%ZZ = getPartitionFunction( II{2}, Ev, Tv );
%6*(rmax+1)^2,
    %Temp comparison
%     h0 = openfig('stretch_comp_xy2_b_10_n_30_maxn_100.fig');
%     [~,img] = h0.Children.Children;
%     Ixy2 = img.CData;
%     limsX = img.XData;
%     limsY = img.YData;
%     
%     xs = (limsX <= 4);
%     ys = (limsY >= 0);
    
    %I versus T and E (II)
    %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    limsE = Ev;
    es = (0 <= Ev & Ev <= 5);
    Tvis = (Tv <= nt*Tc);
   % pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
    pcolor(limsE(es),Tv(Tvis), II{5}(Tvis,es) );
    axis([min(limsE(es)), max(limsE(es)), -inf,inf])
    colorbar
    colormap hot
    shading flat
    m3 = max(max( II{5}(Tvis,es) ));
    caxis([0,m3])
    xlabel('\omega/J');
    ylabel('T/J');
    hold off;
    filename = ['stretch_comp_mcmc_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(floor(1000*s)),'_special'];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    


toc