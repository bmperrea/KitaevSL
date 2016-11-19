%This script runs parallel mcmc runs for different temperatures
startTime = clock;

problem.restart = true; %Whether to start from the last run
problem.plotting = false; %used for real time debugging/diagnostics

problem.rmax = 15; %radius of honeycomb flake in plaquettes
rmax = problem.rmax;
problem.s = 0.04; % * (rmax + 1/2)/(15 + 1/2); %stretch amount
problem.b = 10;   % magnetic grunesien parameter
problem.stat_rate = 200; %How many proposals before recomputing statistics
problem.plot_rate = 1; %how many statistics before plotting
problem.max_batches = 2000;
problem.max_count = 100000; %50,000


% msi instructions
% start a job
% qsub -q small setup.pbs
% cancel a job
% qsig -s SIGINT 2689815.mesabim3.msi.umn.edu


s = problem.s;
b = problem.b;

Trel = [0.001,.01,.1,.25,.35,.4,.45,.5,.55,.575,.6,.625,.65,.675,.7,.725,...
    .75,.775,.8,.825,.85,.875,.9,.925,.95,...
    1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.8,2,5,10,...
    20,30,40,50,60,70,80,100,1000,10000];
Tcomp = [.15,3,29.2,.03,.054,.078]; 

numplaq = 3*problem.rmax*(problem.rmax + 1) + 1;
Tc = .13/log(numplaq); disp(Tc)

Tv = sort( [Tc*Trel, Tcomp] ); 

%Tv = Tv(Tv>0.05)

%Tv = Tv(1:5);
emax = 12.0; bins = 200;
Ev = (1:bins)'*emax/bins; limsE = Ev;

format compact
%close all

%Main Loop
disp(['rmax = ',num2str(rmax)])
disp(numel(Tv))

%Shuffle the indices to help the task distribution be more even
% on the parallel cores
order = randperm(numel(Tv));

outputs = cell(1,numel(Tv));
for j = (1:numel(Tv))
   tic
   outputs{order(j)} = flux_controller_MCMC(problem,Tv(order(j)));
   disp(order(j))
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

%     %I versus T and E (II)
%     hh=figure;%('Position',position);
%     hold on;
%     pcolor(limsE,Tv(Tvis),II{2}(Tvis,:));
%     axis([0, emax, -inf, inf])
%     colorbar
%     shading flat
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
%     pcolor(limsE,Tv(Tvis),(II{3}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     shading flat
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
%     pcolor(limsE,Tv(Tvis),(II{4}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     shading flat
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
%     pcolor(limsE,Tv(Tvis),(II{5}(Tvis,:)));
%     axis([0, emax, -inf, inf])
%     colorbar
%     shading flat
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixy2_mcmc_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
% 
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

    %Save the data for plotting later
    save(['stretch_flux_mcmc1_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(round(1000*s))]...
        ,'Tc','Tv','Tvis','Ev','II','Ier','outputs','Es','ps')
    
stopTime = clock;
disp('Time in hours:')
disp(etime(stopTime,startTime)/3600);