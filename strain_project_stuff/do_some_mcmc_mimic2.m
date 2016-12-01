%This script runs parallel mcmc runs for different temperatures
startTime = clock;

problem.restart = true; %Whether to start from the last run
problem.plotting = true; %used for real time debugging/diagnostics

problem.rmax = 11; %radius of honeycomb flake in plaquettes
problem.s = 0.00; %stretch amount
problem.b = 10;   % magnetic grunesien parameter
problem.stat_rate = 200; %How many proposals before recomputing statistics
problem.plot_rate = 1; %how many statistics before plotting
problem.max_batches = 2000;
problem.max_count = 50000; %50,000

rmax = problem.rmax;
s = problem.s;
b = problem.b;

Tv =[0,.15,3,29.2,100000]; 
%[.15,3,29.2
Tc = .13/log(numplaq); disp(Tc)

%Tv = Tv(1:5);
emax = 12.0; bins = 200;
Ev = (1:bins)'*emax/bins; limsE = Ev;

format compact
close all

%Main Loop

for j = 1:numel(Tv) 
   tic
   outputs{j} = flux_controller_MCMC(problem,Tv(j));
   disp(j)
   toc
end

%Aggregate the data
II = cell(1,5);
Ier = cell(1,5);
for u = 2:5
    II{u} = zeros(numel(Tv),bins);
    Ier{u} = zeros(numel(Tv),bins);
end
Es = zeros(1,numel(Tv));
ps = zeros(1,numel(Tv));
for j = 1:numel(Tv) 
    II{2}(j,:) = outputs{j}.dd.';
    II{3}(j,:) = outputs{j}.Ixx.';
    II{4}(j,:) = outputs{j}.Ixy.';
    II{5}(j,:) = outputs{j}.Ixy2.';
    %
    Ier{2}(j,:) = outputs{j}.dde;
    Ier{3}(j,:) = outputs{j}.Ixxe;
    Ier{4}(j,:) = outputs{j}.Ixye;
    Ier{5}(j,:) = outputs{j}.Ixy2e;
    %
    Es(j) = outputs{j}.En;
    ps(j) = outputs{j}.p;
end


%% Output
savePlots = true;

%Temperatures to plot against (skip the high temps to keep it easy to
        %look at)
Tvis = (Tv <= 3*Tc);
prep_plots;

figure;
plot(Ev,II{4})

%   %I versus T and E (II)
%     hh=figure;%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis),II{2}(Tvis,:));
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_DOS_mcmc_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     
%       %I versus T and E (II)
%     hh=figure;%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis),(II{3}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixx_mcmc_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     
%       %I versus T and E (II)
%     hh=figure;%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis),(II{4}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixy_mcmc_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     
%       %I versus T and E (II)
%     hh=figure;%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis),(II{5}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixy2_mcmc_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     %Energy versus T
%     hh=figure;
%     hold on;
%     plot(Tv(Tvis),Es(Tvis))
%     hold off;
%     filename = ['stretch_En_mcmc_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     %p versus T
%     hh=figure;
%     hold on;
%     plot(Tv(Tvis),ps(Tvis))
%     hold off;
%     filename = ['stretch_ps_mcmc_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
% 
%     %Save the data for plotting later
    save(['stretch_flux_mcmc_mimic2_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)]...
        ,'Tc','Tv','Tvis','Ev','II','Ier','outputs','Es','ps')
    
stopTime = clock;
disp('Time in hours:')
disp(etime(stopTime,startTime)/3600);