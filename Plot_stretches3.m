prep_plots
close all
savePlots = true;

favColor1 = [153,216,201]/255 - .25;
favColor2 = [44,152,105]/255 - .15;
lw = 2;

% Load data
rmax = 7;
b = 10;
s=0;

fileName = ['stretch_flux_mcmc_mimic1_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
load(['../',fileName])
%ZZ = getPartitionFunction( 6*(rmax+1)^2, II{2}, Ev, Tv );

h0 = openfig('../2D_DOS_10^7_Jx_100over100.fig');
D0 = h0.Children.Children.YData;

emax = max(Ev);
dE = emax/numel(Ev);

% Plot stuff
Tcomp = [0.1500 3.0000 29.2000 0.0600 0.1080 0.1560];
Tcinds= [31, 39, 43, 8, 22, 27]; %[31, 39, 43, 8, 22, 27];

set(0,'defaultAxesFontSize',13);
 %DOS versus E (II{2})
 %This plot is to compare to Fig 3(a) of PhysRevB.92.115122
    hh=figure;%('Position',position);
    hold on;
    area(Ev/2,D0,'FaceColor',.82*[1,1,1],'EdgeColor',.82*[1,1,1]);
%     plot(Ev/2,normalate( II{2}(Tcinds(4),:) , 2/dE ),':');
%     plot(Ev/2,normalate( II{2}(Tcinds(5),:) , 2/dE ),'-');
%     plot(Ev/2,normalate( II{2}(Tcinds(6),:) , 2/dE ),'-');
%     plot(Ev/2,normalate( II{2}(numel(Tv),:) , 2/dE ),'-');
    plot(Ev/2, II{2}(2,:) ,':')
    plot(Ev/2, II{2}(3,:) ,'-')
    plot(Ev/2, II{2}(4,:) ,'-')
    plot(Ev/2, II{2}(5,:) ,'-')
axis([0, emax/2, -inf, 0.75])
    caxis auto
    xlabel('\omega');
    ylabel('DOS');
    legend({'T=0', ['T=',num2str(Tv(2))], ['T=',num2str(Tv(3))],...
        ['T=',num2str(Tv(4))],'T=\infty'}, 'Location', 'NorthEast');
    hold off;
    filename = ['stretch_DOS_special_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
rmax = 11;
fileName =['stretch_flux_mcmc_mimic2_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
load(['../',fileName])

%2D_DOS_10^7_Jx_100over100
%2D_Raman_xx_10^7_Jx_100over100

h0 = openfig('../2D_Raman_xx_10^7_Jx_100over100.fig');
Ixx0 = h0.Children.Children.YData;
Ixx000 = Ixx0;

bin = 2:200;
%Comparison to 1602.05277v1 Fig 3(a)
hh=figure;%('Position',position);
    hold on;
    area(Ev(bin), normalate( Ixx0(bin), 1/dE ) ,'FaceColor',.82*[1,1,1],'EdgeColor',.82*[1,1,1]);
    plot(Ev(bin), II{4}(2,bin) ,...
        Ev(bin), II{4}(3,bin )...
        ,Ev(bin), II{4}(4,bin) );
    axis([0, emax, -inf, .24])
 %   caxis auto
    xlabel('\omega');
    ylabel('I_{xx}');
    legend({'T=0', ['T=',num2str(Tv(2))], ['T=',num2str(Tv(3))],...
        ['T=',num2str(Tv(4))]}, 'Location', 'NorthEast');
    hold off;
    filename = ['stretch_Ixx_special_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
    %colorplots 
    
    %Size comparison
    h0 = openfig('size_comp_xy2_b_10_s_4_maxn_50.fig');
    [~,img] = h0.Children.Children;
    Ixy2 = img.CData;
    limsX = img.XData;
    limsY = img.YData;    
    
    save('size_comp_xy2_b_10_s_4_maxn_50','limsX','limsY','Ixy2') 
    
    xs = (limsX <= 4);
    ys = (limsY >= 10);

    
    set(0,'defaultAxesFontSize',18);
    
    %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    pcolor(limsX(xs),limsY(ys),norma(Ixy2(ys,xs),2));
    axis([min(limsX(xs)), max(limsX(xs)), min(limsY(ys)), max(limsY(ys))])
    colorbar
    colormap hot
    shading flat
    caxis auto
    xlabel('\omega');
    ylabel('r');
    hold off;
    filename = ['size_comp_xy2_b_10_s_4_maxn_50_special'];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
    
    %Stretch comparison
    h0 = openfig('stretch_comp_xy2_b_10_n_30_maxn_100.fig');
    [~,img] = h0.Children.Children;
    Ixy2 = img.CData;
    limsX = img.XData;
    limsY = img.YData;
    
    save('stretch_comp_xy2_b_10_n_30_maxn_100','limsX','limsY','Ixy2')   
        
    xs = (limsX <= 4);
    ys = (limsY >= 0);
    
    %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    pcolor(limsX(xs),limsY(ys),norma(Ixy2(ys,xs),2));
    axis([min(limsX(xs)), max(limsX(xs)), min(limsY(ys)), max(limsY(ys))])
    colorbar
    colormap hot
    shading flat
    caxis auto
    xlabel('\omega');
    ylabel('s');
    hold off;
    filename = ['stretch_comp_xy2_b_10_n_30_maxn_100_special'];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
    
    
rmax = 18;
s = .04;
fileName = ['stretch_flux2_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
load(['../',fileName])
disp(Ev)
ZZ = getPartitionFunction( II{2}, Ev, Tv );
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
    Tvis = (Tv <= 2*Tc);
   % pcolor(limsE(es),Tv(Tvis),norma(II{5}(Tvis,es),2));
    pcolor(limsE(es),Tv(Tvis), II{5}(Tvis,es) ./ repmat(ZZ(Tvis),1,sum(es)) );
    axis([min(limsE(es)), max(limsE(es)), -inf,inf])
    colorbar
    colormap hot
    shading flat
    m3 = max(max( II{5}(Tvis,es) ./ repmat(ZZ(Tvis),1,sum(es)) ));
    caxis([0,m3])
    xlabel('\omega');
    ylabel('T');
    hold off;
    filename = ['stretch_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s),'_special'];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
     %Energy versus p
     
    nup = sum(pvs);
    pvals = zeros(1,nup);
    Evalys = pvals;
    ii = 1;
    for ind = 1:numel(ps)
        pind = ps(ind);
        for tt = 1:pvs(pind)
            pvals(ii) = pind-1;
            Evalys(ii)= Endist{pind}(tt);
            ii = ii+1;
        end
    end
     
    hh=figure;
    hold on;
    scatter(pvals/numel(pvs),(Evalys-Evalys(1))/numel(pvs),1,'filled')
    xlabel('p/N_P');
    ylabel('Energy per plaq.');
    hold off;
    filename = ['stretch_En_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_s_',num2str(round(1000*s))];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Basic plots
    
   % rmax = 50;
   % s = .04;
    h0 = openfig('stretch_I2_rmax_50_b_10_c_40_p_0');
    [~,dat] = h0.Children.Children;
    Ixx = dat(11).YData;
    Ev = dat(11).XData;
    
    h0 = openfig('stretch_DOS_rmax_50_b_10_c_40_p_0');
    Ddd = h0.Children.Children(10).YData;
    
    h0 = openfig('stretch_I3_rmax_50_b_10_c_40_p_0');
    [~,dat] = h0.Children.Children;
    Ixy2 = dat(10).YData;
    
    
    h0 = openfig('stretch_I2_rmax_50_b_10_c_0_p_0');
    [~,dat] = h0.Children.Children;
    Ixx0 = dat(11).YData;
    Ev = dat(11).XData;
    
    h0 = openfig('stretch_DOS_rmax_50_b_10_c_0_p_0');
    Ddd0 = h0.Children.Children(10).YData;
    
    h0 = openfig('stretch_I3_rmax_50_b_10_c_0_p_0');
    [~,dat] = h0.Children.Children;
    Ixy20 = dat(10).YData;
    
    
    s=0.04; rmax = 50;
    c = s/(sqrt(3)*rmax+sqrt(3)/2);
    Bfield = -8*c*b;
    E0 = 3*sqrt(abs(Bfield));
    peaks = 7;

    %Plot DOS
    hh=figure;
    hold on;
    plot(Ev/2,Ddd0,':','Color',favColor1,'LineWidth',lw)
    plot(Ev/2,Ddd,'Color',favColor2,'LineWidth',lw);
    xs = (Ev/2 <= 2);
    ys = ones(size(Ddd));
 %   xlim([0,6.45])
    %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
    %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
    xlabel('\omega');
    ylabel('DOS');
    legend({'s=0','s=0.04'}, 'Location', 'NorthEast');
    axis([min(Ev(xs)/2), max(Ev(xs)/2), 0, max(Ddd)])
    %set(gca,'XTick',-3:3); 
    %set(gca,'YTick',2*(0:5));
    %hold off;
    yl = ylim;
    for n = 0:peaks
        %[Xi, ~] = ds2nfu(sqrt(n)*E0, 0);
        %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')    
        line(sqrt(n)*[E0,E0],yl,'Color',[.7,.7,.7])
    end
%    plot(Ev/2,Ddd,Ev/2,Ddd0);
    hold off;
    filename = ['stretch_DOS_special_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);

    %Plot the ramans
    hh=figure;%('Position',position);
    hold on;
    plot(Ev,Ixy20,':','Color',favColor1,'LineWidth',lw)
    plot(Ev,Ixy2,'Color',favColor2,'LineWidth',lw);
    xs = (Ev <= 4);
    %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
    %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
    xlabel('\omega');
    ylabel('I_{[xy]}');
    legend({'s=0','s=0.04'}, 'Location', 'NorthEast');
    axis([min(Ev(xs)), max(Ev(xs)), 0, max(Ixy2)])
    yl = ylim;
    for n = 0:peaks
        %[Xi, ~] = ds2nfu(2*sqrt(n)*E0, 0);
        %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1') 
        line(2*sqrt(n)*[E0,E0],yl,'Color',[.7,.7,.7])
    end
 %   plot(Ev,Ixy2,Ev,Ixy20);
    hold off;
    filename = ['stretch_I3_special_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);

    hh=figure;%('Position',position);
    hold on;
    plot(Ev,Ixx0,':','Color',favColor1,'LineWidth',lw);
    plot(Ev,Ixx,'Color',favColor2,'LineWidth',lw);
    xs = (Ev <= 4);
    %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
    %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
    xlabel('\omega');
    ylabel('I_{(xx)}');
 %   legend({'I_{xx}'}, 'Location', 'NorthEast');
    legend({'s=0','s=0.04'}, 'Location', 'NorthEast');
    axis([min(Ev(xs)), max(Ev(xs)), 0, max(Ixx)])
    yl = ylim;
    for n = 0:peaks-1
        %[Xi, ~] = ds2nfu((sqrt(n)+sqrt(n+1))*E0, 0);
        %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')  
        line((sqrt(n) + sqrt(n+1))*[E0,E0],yl,'Color',[.7,.7,.7])
    end
    %plot(Ev,Ixx,Ev,Ixx0);
    hold off;
    filename = ['stretch_I2_special_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);

    
    %Full spectrum comparison
     hh=figure;%('Position',position);
    hold on;
    plot(Ev*12/12.5,Ixx000 / 2,'Color',[0,0,0],'LineWidth',lw);
 %   plot(Ev,Ixx0,'Color',favColor1,'LineWidth',lw);
    plot(Ev,Ixx,'Color',favColor2,'LineWidth',lw);
    xs = true(size(Ev));
    %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
    %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
    xlabel('\omega');
    ylabel('I_{(xx)}');
 %   legend({'I_{xx}'}, 'Location', 'NorthEast');
    legend({'s=0 (Bulk)','s=0.04'}, 'Location', 'NorthEast');
    axis([min(Ev(xs)), 12, 0, inf])
 %   yl = ylim;
  %  for n = 0:peaks-1
        %[Xi, ~] = ds2nfu((sqrt(n)+sqrt(n+1))*E0, 0);
        %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')  
 %       line((sqrt(n) + sqrt(n+1))*[E0,E0],yl,'Color',[.7,.7,.7])
 %   end
    %plot(Ev,Ixx,Ev,Ixx0);
    hold off;
    filename = ['stretch_I2_special_full_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);
   
    
        %Full spectrum comparison
     hh=figure;%('Position',position);
    hold on;
    plot(Ev*12/12.5,Ixx000 / 2,'Color',[0,0,0],'LineWidth',lw);
    plot(Ev,Ixx0,':','Color',favColor1,'LineWidth',lw);
    plot(Ev,Ixx,'Color',favColor2,'LineWidth',lw);
    xs = true(size(Ev));
    %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
    %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
    xlabel('\omega');
    ylabel('I_{(xx)}');
 %   legend({'I_{xx}'}, 'Location', 'NorthEast');
    legend({'Bulk','s=0','s=0.04'}, 'Location', 'NorthEast');
    axis([min(Ev(xs)), 12, 0, max(Ixx)*.8])
 %   yl = ylim;
  %  for n = 0:peaks-1
        %[Xi, ~] = ds2nfu((sqrt(n)+sqrt(n+1))*E0, 0);
        %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')  
 %       line((sqrt(n) + sqrt(n+1))*[E0,E0],yl,'Color',[.7,.7,.7])
 %   end
    %plot(Ev,Ixx,Ev,Ixx0);
    hold off;
    filename = ['stretch_I2_special_full2_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);
   