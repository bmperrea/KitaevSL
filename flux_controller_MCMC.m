prep_plots;

%These initial parameters define the problem/task for the code
rmax = 6;
s = 0.04;
b = 10;
flag = -1; %(don't plot)
nmax = 5;
T = 0.015;

batch_size = 10;
plot_rate = 10;

prob_rand = 1/4;
prob_shuffle = 1/4;
prob_move  = 1/4;

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
p0 = 0; %Could use logic to find a good start! Otherwise this uses 
            %reverse annealing
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

%Initialize a bunch of storage variables.
count = 0;
decision = 0;
happy = false;

%Declare plots
    h1= figure;
    title('1/ max mean err')
    xlabel('cycles');
    ylabel('Inverse Max (over observable) mean (over E) Error');   

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
    
    if decision == 1 || decision == 4
        % Convert plaquette index data to polar index coordinates
        [ru,ku,mu] = get_polar(plaq);

        % Do the actual flipping to a temporary Hamiltonian
        [H1,Rxx1,Rxy1,RxyA1,RxyB1] = flip_plaquette(ru,ku,mu,rmax,...
                                                    H,Rxx,Rxy,RxyA,RxyB);
        pbool(plaq) = ~pbool(plaq);
    end

    %if moving a plaq, also do the removal
    if decision == 3
        [ru2,ku2,mu2] = get_polar(plaq2);
        [H1,Rxx1,Rxy1,RxyA1,RxyB1] = flip_plaquette(ru2,ku2,mu2,rmax,...
                                        H1,Rxx1,Rxy1,RxyA1,RxyB1);
        pbool(plaq2) = ~pbool(plaq2);       
    end
    
    if decision == 2
        %Take a new random set of fluxes with same number.
        p = sum(pbool);
        [H1,Rxx1,Rxy1,RxyA1,RxyB1,pbool] = stretch_2D_makeH(rmax,b,s,p);        
    end
    
    
%% Test

    En1 = -sum(svd(H1)); %This is time consuming.
    dEn = En1 - En0;
    
    %Decide whether to accept the proposal state
    accept = false;
    if dEn < 0
        accept = true;
    else 
        prob = exp(-dEn/T);
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
        %actually compute the observables of interest
        [~,~,dd,~,Ixx,Ixy,Ixy2] = ...
            dos1_loop(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax);
    else
        %Just add the same state again in the update! (Do nothing)
    end
    
%% Update

    if mod(count,batch_size) == 0

        %add this data to the Markov Chains (there is only one state chain, 
            % but several sets of value chains).

        dds   = [ dds  ; reshape(dd  ,1,bins) ];
        Ixxs  = [ Ixxs ; reshape(Ixx ,1,bins) ];
        Ixys  = [ Ixys ; reshape(Ixy ,1,bins) ];
        Ixy2s = [ Ixy2s; reshape(Ixy2,1,bins) ];
        Ens =   [ Ens  ; En ];
        decisions = [decisions;decision];
        happys = [happys; happy];

        %Do some error analysis and update the error chains.

        [ddt,ddt1]     = geyer_icse(dds);
        [Ixxt,Ixxt1]   = geyer_icse(Ixxs);
        [Ixyt,Ixyt1]   = geyer_icse(Ixys);
        [Ixy2t,Ixy2t1] = geyer_icse(Ixy2s);
        [Ent,Ent1]     = geyer_icse(Ens);

        ddt2   = geyer_imse(dds);
        Ixxt2  = geyer_imse(Ixxs);
        Ixyt2  = geyer_imse(Ixys);
        Ixy2t2 = geyer_imse(Ixy2s);
        Ent2   = geyer_imse(Ens);  

        ddst   = [ ddst  ; reshape(ddt  ,1,bins) ];
        Ixxst  = [ Ixxst ; reshape(Ixxt ,1,bins) ];
        Ixyst  = [ Ixyst ; reshape(Ixyt ,1,bins) ];
        Ixy2st = [ Ixy2st; reshape(Ixy2t,1,bins) ];
        Enst =   [ Enst  ; Ent ];

        ddst1   = [ ddst1  ; reshape(ddt1  ,1,bins) ];
        Ixxst1  = [ Ixxst1 ; reshape(Ixxt1 ,1,bins) ];
        Ixyst1  = [ Ixyst1 ; reshape(Ixyt1 ,1,bins) ];
        Ixy2st1 = [ Ixy2st1; reshape(Ixy2t1,1,bins) ];
        Enst1 =   [ Enst1  ; Ent1 ];

        ddst2   = [ ddst2  ; reshape(ddt2  ,1,bins) ];
        Ixxst2  = [ Ixxst2 ; reshape(Ixxt2 ,1,bins) ];
        Ixyst2  = [ Ixyst2 ; reshape(Ixyt2 ,1,bins) ];
        Ixy2st2 = [ Ixy2st2; reshape(Ixy2t2,1,bins) ];
        Enst2 =   [ Enst2  ; Ent2 ];
        
        % Test happiness
        ratios = [mean(sqrt(ddst(end,:)))/mean(dds),...
                    mean(sqrt(Ixxst(end,:)))/mean(Ixxs),...
                    mean(sqrt(Ixyst(end,:)))/mean(Ixys),...
                    mean(sqrt(Ixy2st(end,:)))/mean(Ixy2s)];
        if max(ratios) < lowenough
            happy = true;
        end

        %display some stuff
        disp(count)
        disp(max(ratios))    

        ratios1 = [mean(sqrt(ddst1(end,:)))./mean(dds),...
                    mean(sqrt(Ixxst1(end,:)))./mean(Ixxs),...
                    mean(sqrt(Ixyst1(end,:)))./mean(Ixys),...
                    mean(sqrt(Ixy2st1(end,:)))./mean(Ixy2s)];
        ratios2 = [mean(sqrt(ddst2(end,:)))./mean(dds),...
            mean(sqrt(Ixxst2(end,:)))./mean(Ixxs),...
            mean(sqrt(Ixyst2(end,:)))./mean(Ixys),...
            mean(sqrt(Ixy2st2(end,:)))./mean(Ixy2s)];
        
        if mod(count,batch_size*plot_rate) == 0
            %Plot some error diagnostics
            hh = figure(h1);
            plot();     
            plot();
            plot();
        end     
        
    end    
    disp(accept)
    
end


%% Post-process

% Throw out annealing stage to get a better sample mean
    % This can be the initial sequence up to the first time the energy is
    % within "1 std" of the mean. Here std is ~sqrt(N)*sqrt(var)
Estd  = sqrt(Ent) * sqrt( batch_size*numel(Ent) );
reasonables = find(Estd < mean(Ens) + Estd);
if reasonables(1) > numel(Ent)/2
    warning('It took forever to anneal?')
end
start_err = floor(reasonables(1)/batch_size);
start = start_err * batch_size + 1;
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


