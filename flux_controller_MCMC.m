% Plot key
%
% 1/error        energy           flux sector p
% decisions      acorr time       raw sample var

restart = false;
plotting= true;

%These initial parameters define the problem/task for the code
rmax = 15;
s = 0.04;
b = 10;
flag = -1; %(don't plot)
T = 0.015;


stat_rate = 100;
plot_rate = 1;

prob_rand = 1/4;
prob_shuffle = 1/16; %Turns out shuffle proposals are slow.
prob_move  = 1/4;
% prob_NNN = 7/16

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes


startTime = clock;
lowenough = 1/(65); %65 appears to be the default number of bins in the default colormap in my version of matlab

%Whether to plot 
plotting = true;

%save the command window to text file
fileN = ['stretch_mcmc_diary_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
diary(fileN)
diary on

%number of plaquettes
numplaq = 1 + 3*rmax*(rmax+1);

%some other choices
knum = 5;
emax = 12.0; bins = 200;
Ev = (1:bins)'*emax/bins;
emax = max(Ev);

%Get the base Hamiltonian
H = stretch_2D_makeH(rmax,b,s,0);
En00 = -sum(svd(H)); %Tracks the ground state energy

% Set a base gauge 
try %start from previous run
    load(['stretch_mcmc_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s),'_T_'...
        num2str(round(T*1000))])
catch 
    restart = true;
end

if restart
       
    p0 = 2; %Could use logic to find a good start! Otherwise this uses 
            %reverse annealing
    disp(['no initial data, starting from p0 = ',num2str(p0)])
    
    [H,Rxx,Rxy,RxyA,RxyB,pbool,En0] = stretch_2D_makeH(rmax,b,s,p0);
    %En0 Tracks the current state energy
    En = En0 - En00; %Tracks the current state excitation energy 

    %Take the first data
    [~,~,dd,~,Ixx,Ixy,Ixy2] = ...
                dos1_loop(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax);

    dds   = reshape(dd  ,1,bins);
    Ixxs  = reshape(Ixx ,1,bins);
    Ixys  = reshape(Ixy ,1,bins);
    Ixy2s = reshape(Ixy2,1,bins);
    Ens = En;
    decisions = 0;
    happys = false;
    changed = false;
    pss = p0;

    ddst   = dds;
    Ixxst  = Ixxs;
    Ixyst  = Ixys;
    Ixy2st = Ixy2s;
    Enst = Ens;

    raw_vars = [];
    errors   = [];
    acorrtimes = [];

    %Initialize a bunch of storage variables.
    count = 0;
    decision = 0;
    happy = false;
end



%Declare plots
%    h1= figure;
%    title('1/ max mean err')
%    xlabel('cycles');
%    ylabel('Inverse Max (over observable) mean (over E) Error');   
%    hold on

if plotting    
    prep_plots;
    close all;
    h0 = figure('units','normalized','outerposition',[0 0 1 1]);
end

%h1 = figure;

%h2 = figure;

% MCMC main loop
while ~happy

    count = count + 1;
    
%% Proposal


    ra = rand;
    if sum(pbool) == 0 || ra < prob_rand
        
        decision = 1;
        % flip a random plaquette
        plaq = randi(numplaq);
        
    elseif ra < prob_rand + prob_shuffle
                
        decision = 2;
        %Shuffle the fluxes randomly
        
    elseif 1-ra < prob_move       
        
        decision = 3;
        
        % move a random plaquette to a random plaquette
        fluxes = find(pbool);
        removal = randi(numel(fluxes));
        plaq2 = fluxes(removal);
        
        nfluxes = find(~pbool);
        insertion = randi(numel(nfluxes));
        plaq = nfluxes(insertion);
        
    else % flip a 2nd neighbor plaquette
        
        decision = 4;
        %find the second neighbors of the fluxes
        fluxes = find(pbool);
        neighbors = fluxes;
        for ind = 1:numel(fluxes)
            flux = fluxes(ind);   
            %Convert plaquette number to polar coordinate data
            [r,k,m] = get_polar(flux);            
            %get the neighbors
            nmax = 2;
            ps = get_NNNNs(r,k,m,nmax);
            ps = ps(ps<=numplaq); %remove the stuff outside the lattice
            %add them in, removing duplicates
            neighbors = unique( [neighbors,ps] );
        end
        
        %pick one
        ind = randi(numel(neighbors));
        plaq = neighbors(ind);
            
    end
    
    if decision ~= 2
        % Convert plaquette index data to polar index coordinates
        [ru,ku,mu] = get_polar(plaq);

        % Do the actual flipping to a temporary Hamiltonian
        [H1,Rxx1,Rxy1,RxyA1,RxyB1] = flip_plaquette(ru,ku,mu,rmax,...
                                                    H,Rxx,Rxy,RxyA,RxyB);
        pbool1 = pbool;
        pbool1(plaq) = ~pbool1(plaq);
    end

    %if moving a plaq, also do the removal
    if decision == 3
        [ru2,ku2,mu2] = get_polar(plaq2);
        [H1,Rxx1,Rxy1,RxyA1,RxyB1] = flip_plaquette(ru2,ku2,mu2,rmax,...
                                        H1,Rxx1,Rxy1,RxyA1,RxyB1);
        pbool1(plaq2) = ~pbool1(plaq2);
    end
    
    if decision == 2
        %Take a new random set of fluxes with same number.
        p = sum(pbool);
        [H1,Rxx1,Rxy1,RxyA1,RxyB1,pbool1] = stretch_2D_makeH(rmax,b,s,p);         
    end
    
    
%% Test

    En1 = -sum(svd(H1))-En00; %This is time consuming.
    dEn = En1 - En;
    
    %Decide whether to accept the proposal state
    accept = false;
    
    if decision == 4 %correct for different entropies of 2nd flux configs
                     %This makes the proposal distribution reversible
        num = numel(ps);
        
        % find number of new neighbors
        numel(ps);
        
        fluxes1 = find(pbool1);
        neighbors1 = fluxes1;
        for ind = 1:numel(fluxes1)
            flux = fluxes1(ind);   
            %Convert plaquette number to polar coordinate data
            [r,k,m] = get_polar(flux);            
            %get the neighbors
            nmax = 2;
            ps1 = get_NNNNs(r,k,m,nmax);
            ps1 = ps1(ps1<=numplaq); %remove the stuff outside the lattice
            %add them in, removing duplicates
            neighbors1 = unique( [neighbors1,ps1] );
        end
        
        num1 = numel(ps1); 
        
        dEn = dEn - T*log(num/num1); 
    end
        
    if dEn < 0
        accept = true;
    else 
        prob = exp(-dEn/T);
       % disp(prob)
        ran  = rand;
        if ran < prob
            accept = true; 
        end
    end
    
    if accept
        H    = H1;
        Rxx  = Rxx1;
        Rxy  = Rxy1;
        RxyA = RxyA1;
        RxyB = RxyB1;
        pbool = pbool1;
        En = En1;
        %actually compute the observables of interest
        [~,~,dd,~,Ixx,Ixy,Ixy2] = ...
            dos1_loop(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax);
        
        changed = true; %Whether it has changed from the start at all;
                
       % disp(sum(pbool))
    else
        decision = 0;
        %Just add the same state again in the update! (Do nothing)
    end
    
%% Update

    p = sum(pbool);

    %add this data to the Markov Chains (there is only one state chain, 
        % but several sets of value chains).

    dds   = [ dds  ; reshape(dd  ,1,bins) ];
    Ixxs  = [ Ixxs ; reshape(Ixx ,1,bins) ];
    Ixys  = [ Ixys ; reshape(Ixy ,1,bins) ];
    Ixy2s = [ Ixy2s; reshape(Ixy2,1,bins) ];
    Ens =   [ Ens  ; En ];
    pss  =  [ pss  ; p  ];
    decisions = [decisions;decision];
    happys = [happys; happy];

    if mod(count,stat_rate) == 0 && changed

        %Do some error analysis and update the error chains.

        ddt    = initseq_vec(dds.');
        Ixxt   = initseq_vec(cast(Ixxs.','double'));
        Ixyt  = initseq_vec(cast(Ixys.','double'));
        Ixy2t = initseq_vec(cast(Ixy2s.','double'));
        Ent   = initseq_vec(Ens.');  

        ddst   = [ ddst  ; reshape(ddt  ,1,bins) ];
        Ixxst  = [ Ixxst ; reshape(Ixxt ,1,bins) ];
        Ixyst  = [ Ixyst ; reshape(Ixyt ,1,bins) ];
        Ixy2st = [ Ixy2st; reshape(Ixy2t,1,bins) ];
        Enst   = [ Enst  ; Ent ];
        
        % Test happiness
        stds = [sqrt(mean(ddst.'));...
                    sqrt(mean(Ixxst.'));...
                    sqrt(mean(Ixyst.'));...
                    sqrt(mean(Ixy2st.'))].';
                
        means = [max(mean(dds)),max(mean(Ixxs.')),max(mean(Ixys.')),...
                    max(mean(Ixy2s.'))];
                
        error = ( stds(end,:)./means ) ./ sqrt(count);
                
        errors = [errors ; error];
        
        raw_var = [mean(var(dds)),mean(var(Ixxs)),...
                    mean(var(Ixys)),mean(var(Ixy2s))];
                
        acorrtime = stds(end,:).^2./raw_var;
        acorrtimes = [acorrtimes;acorrtime];
                
        raw_vars = [raw_vars ; raw_var./means];
        
        if max(error) < lowenough
            happy = true;
        end

        %display some stuff
       % disp(sum(pbool))
       % disp(max(errors))    
         
        
        %Save the data everytime we compute error
        try
            save(['stretch_mcmc_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s),'_T_'...
                num2str(round(T*1000))]...
            ,'dds','Ixxs','Ixys','Ixy2s','Ens','pss','decisions','happys',...
            'ddst','Ixxst','Ixyst','Ixy2st','Enst','ratios','errors',...
            'error','stds','acorrtimes','raw_vars','happy',...
            'pbool','count','H','Rxx','Rxy','RxyA','RxyB')
        catch
            disp('couldnt save')
        end
        
        
        if mod(count,stat_rate*plot_rate) == 0 && plotting
            %Plot some error diagnostics
            hh = figure(h0); %hold on;
            subplot(2,3,1)
            title('1/ max mean err')
           % xlabel('cycles');
          %  ylabel('Inverse Max (over observable) mean (over E) Error');
            %hold on
            %plot([0,count/stat_rate],[65,65]);  
            plot(1./errors);   
            refline(0,1/lowenough);
            %hold off
            
            %Plot some other things
        %    hh = figure(h1); %hold on;
        subplot(2,3,2)
            title('Energy')
        %    xlabel('cycles');
        %    ylabel('E/J');   
            plot(Ens);   
            
        %    hh = figure(h2); %hold on;
        subplot(2,3,3)
            title('flux sector')
            xlabel('cycles');
            ylabel('p');   
            plot(pss);  
            
          %  hh = figure(h3);
          subplot(2,3,4)
         %   title('sample variances')
        %    xlabel('cycles');
        %    ylabel('variance');   
         %   plot(stds);       
            
        %    hh = figure(h1);
        %     subplot(3,1,1)
%             title('decisions')
%             xlabel('cycles');
%             ylabel('decision');   
%             plot(decisions);     
          
            histogram(decisions,'Normalization','probability');
            
            subplot(2,3,5)
           % hh = figure(h2);
            plot(acorrtimes)
            
            subplot(2,3,6)
            plot(raw_vars);
            
            pause(.001) %Make sure the plots show up.
        end    
        
        his = histogram(decisions,'Normalization','probability');
            
        %if we are in a high flux regime 
        %increase probability of shuffle.
%         if p > 0.05*numplaq && (his.Values(3))/(1-his.Values(1)) > 0.25*0.1
%             prob_shuffle = 1/2;
%             prob_move = 1/8;
%             %Then prob NNN = 1/8;
%         end
        
    end    
    
   % disp(accept)
    
end


%% Post-process

% Throw out annealing stage to get a better sample mean
    % This can be the initial sequence up to the first time the energy is
    % within "1 std" of the mean. Here std is ~sqrt(N)*sqrt(var)
Estd  = sqrt(En) * sqrt( stat_rate*numel(En) );
reasonables = find(Estd < mean(Ens) + Estd);
if reasonables(1) > count/2
    warning('It took forever to anneal?')
end
start_err = floor(reasonables(1)/stat_rate);
start = start_err * stat_rate + 1;
start_err = start_err + 1;

dds2   = dds(start:end,:);
Ixxs2  = Ixxs(start:end,:);
Ixys2  = Ixys(start:end,:) ;
Ixy2s2 = Ixy2s(start:end,:);
Ens2 =   Ens(start:end,:)  ;

% compute means
dd = mean(dds2);
Ixx = mean(Ixxs2);
Ixy = mean(Ixys2);
Ixy2 = mean(Ixy2s2);
En   = mean(Ens2);

% Redo error analysis? (Nope)
dde = ddst(end,:);
Ixxe = Ixxst(end,:);
Ixye = Ixxst(end,:);
Ixy2e= Ixy2st(end,:);

if plotting
    % plot some data
    hh = figure;
    plot(Ev,dd);
    title('dos')
    xlabel('\omega/J');
    ylabel('DOS');

    hh = figure;
    plot(Ev,Ixx,Ev,Ixy,Ev,Ixy2);
    title('Raman Intensity')
    xlabel('\omega/J');
    ylabel('I');
end


