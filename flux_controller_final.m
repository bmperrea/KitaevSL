startTime = clock;

% This code controls a very specific calculation of the Raman intensities
% defined in dos1.m by sampling flux sectors set in stretch_flux6.m
% This controller calls stretch_flux6 and computes a running error for each
% temperature in the calculation, continuing only with temperatures that
% have yet to acheive the desired error defined as:
lowenough = 1/(65); %65 appears to be the default number of bins in the default colormap in my version of matlab

%These initial parameters define the problem/task for the code
    %number of honeycomb layers
rmax = 6; %should be between 0 and 18 (out of bound errors occur above 18)
    %magnetic response to stretch
b = 10;
    %max percent stretch
s = 0;

%Whether to plot 
plotting = true;

%save the command window to text file
fileN = ['stretch_flux_diary_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
diary(fileN)
diary on

    %number of plaquettes
numplaq = 1 + 3*rmax*(rmax+1);
    %Temperatures to sample relative to Tc
% Trel = [0.001,.01,.04,.1,.2,.4,.6,.65,.7,.725,.75,.775,.8,.825,.85,.875,.9,.925,.95,.975,1.00,...
%     1.025,1.05,1.075,1.10,1.125,1.15,1.175,1.20,1.225,1.25,1.275,1.3,1.35,1.4,1.6,1.8,2,5,10,...
%     20,30,35,37,39,41,43,45,47,49,51,55,60,80,100,1000,10000];
Trel = [0.001,.01,.1,.25,.35,.4,.45,.5,.55,.575,.6,.625,.65,.675,.7,.725,...
    .75,.775,.8,.825,.85,.875,.9,.925,.95,...
    1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.8,2,5,10,...
    20,30,40,50,60,70,80,100,1000,10000];
%The second transition sould be around 56 to 57 Tc
%temperatures to compare with other papers in units of J
Tcomp = [.15,3,29.2,.03,.054,.078]; 
% indices of these: 31, 38, 43, 8, 22, 27
%First four are relative to 20^2 plaquettes and the last two to 12^2
%That corresponds here to roughly rmax = 11 and rmax = 6

%These initial parameters are tuned for optimal running
    %The minimum number of samples for a given flux sector
rp = 5;
    %The number of points to take in the initial E vs p estimate
pnum = 31; 
    %The number of rp samples to distribute in the first step
pnum2 = 60;
    %The number of rp samples to distribute in later step.
pnum3 = 12;

%Some quick calculations
    %expected Tc
Tc = .13/log(numplaq); disp(Tc)
    %Temperatures to sample
Tv = sort( [Tc*Trel, Tcomp] ); 

disp(Tv)
    
Tvis = (Tv <= 3*Tc);


% Some more initialization

    %Get rid of extra new lines between the outputs
format compact

    %Initialize pool of parallel 'workers' (and track how long it takes to
                                            %do so)
tic
    %If you don't have the Parallel Toolbox you could set parNum to 1 
        %parNum = 1;
    %and comment the following out:
    myCluster = parcluster('local');
    parNum = myCluster.NumWorkers;
    disp(['Working with ',num2str(parNum),' workers'])
    %You may also have to change 'parfor' to 'for' below.
toc

%These initial variables are set in stretch_flux6.m and must be identical
                                                                %here
knum = 5;
emax = 12.0; bins = 200;
Ev = (1:bins)'*emax/bins;
emax = max(Ev);

    %set up font size and shape of plots
prep_plots;
    %Get an initial estimate of the function energy versus p
    % and in doing so an estimate of what flux sectors to probe initially
initialize_pv2;
    %Store the ground state energy, from which to compute Boltzmann factors
    % [If we aren't careful these can easily get beyond numerical
        % precision, which is < 10^(400). ]
GSen = min(Es);
    %make sure we hit the low flux sectors the first time.
pv(1:20) = rp;
    %Set the total pvs to pv initially
pvs = pv;
    %Initialize some variables
        %mean Error per variable (averaged over \omega and T)
meanErrs = [];
        %mean Error per variable (averaged over \omega and max over T)
maxMeanErrs = [];
        %max change per variable relative to the previous step
prevErrs = [];

pf = pf1.';

%prevErrFavs = [];
%meanErrFavs = [];
    
    %initialize the total distribution covered for a given temperature   
pvTemp = cell(1,numel(Tv));

%Now we compute the first step
    % set number of samples
pnum = pnum2;

    % Manually distrubute the problem over parallel workers
    ps0 = find(pv~=0);
    q = numel(ps0);
    list = 1:q;
    parN = min(parNum,q);
    thisList = cell(1,parN);
    pvlist = cell(1,parN);
    for parInd = 1:parN
        pvlist{parInd} = zeros(size(pv));
        thisList{parInd} = ps0( mod(list,parN) == (parInd-1) );
        pvlist{parInd}( thisList{parInd} ) = pv( thisList{parInd} );
    end    
        
    %% The big step
    parfor parInd = 1:parN
        %try
        %disp( pvlist{parInd} )
        [~,IttList{parInd},ItsList{parInd},bfsList{parInd},bfserrList{parInd},EndistList{parInd}] ...
            = stretch_flux6(rmax,b,s,pvlist{parInd},Tv,pf1,GSen); 
        %catch
        %   disp(wtf) 
        %end
    end
    
    %Unpack the distributed results from parallel evaluation
    Itt = IttList{1};
    Its = ItsList{1};
    bfs = bfsList{1};
    bfserr = bfsList{1};
    Endist = EndistList{1};
    pv2 = pvlist{1}; %I compute this just for bug checking
    for parInd = 2:parN
        for kk=1:knum
            Itt( thisList{parInd},kk ) = IttList{parInd}( thisList{parInd},kk );
            Its( thisList{parInd},kk ) = ItsList{parInd}( thisList{parInd},kk );
        end
        bfs( :,thisList{parInd}) = bfsList{parInd}( :,thisList{parInd} );
        bfserr( :,thisList{parInd} ) = bfserrList{parInd}( :,thisList{parInd} );
        Endist( thisList{parInd} ) = EndistList{parInd}( thisList{parInd} );
        pv2( thisList{parInd} ) = pvlist{parInd}( thisList{parInd} );
    end  
    if ~all(pv2 == pv)
        error('screwed up parallel kernel distribution')
    end

pnum = pnum3;

%Set figure handles
close all;

% h1 = figure;
%     title('Raman intensity in I_{[xy]}')
%     xlabel('\omega/J');
%     ylabel('T/J');
% h2 = figure;
%     title('Total Error')
%     xlabel('\omega/J');
%     ylabel('T/J');
% h3 = figure;
%     title('Statistical error')
%     xlabel('\omega/J');
%     ylabel('T/J');
% h4 = figure;
%     title('Integral Error')
%     xlabel('\omega/J');
%     ylabel('T/J');
% h5 = figure;
%     title('error versus p')
%     xlabel('p');
%     ylabel('\delta I');
% h6 = figure;
%     title('statistical err')
%     xlabel('p');
%     ylabel('\delta I_{stat}');
% h7 = figure;
%     title('integration err')
%     xlabel('p');
%     ylabel('\delta I_{int}');
h8 = figure;
    title('scatter of energy values')
    xlabel('p');
    ylabel('E_0/J');
h9 = figure;
    title('1/(mean error)')
    xlabel('cycles');
    ylabel('Inverse mean Error');
h10= figure;
    title('1/ max mean(this one - last one)')
    xlabel('cycles');
    ylabel('Inverse MaxMean difference from previoius');
h11= figure;
    title('1/ max mean err')
    xlabel('cycles');
    ylabel('Inverse Max (over T) mean (over E) Error');
h12= figure;
h13 = figure;
h14=figure;


%MegaLoop
IsayGo = true;
meanVal = 1000000*ones(1,5); %A high number just to initialeze...
 
%Initialize some variables
II = cell(1,knum);
Ier = cell(1,knum); %Ina = cell(1,knum);
Ier1 = Ier; Ier2 = Ier;
Iter = cell(numel(pv),knum);

err1 = cell(1,knum); err2 = err1;
errp1 = cell(1,knum);
errp2 = err1; errp1m = err1; errp2m = err1; errp = err1; errpm = err1;
errI1 = err1; errI2 = err1; errI1m = err1; errI2m = err1;
errI1T = cell(1,5);
errI2T = errI1T;
for kk=1:knum
    %add them in the same loop
    Ier{kk} = 0*Itt{1,kk}; 
    Ier1{kk} = Ier{kk}; Ier2{kk} = Ier{kk};
    II{kk}  = 0*Itt{1,kk};

    err1{kk} = zeros(numel(pv),numel(Tv));
    err2{kk} = err1{kk};

    errp1{kk} = zeros(1,numel(pv));    
    errp2{kk} = zeros(1,numel(pv));

    errp1m{kk} = zeros(1,numel(pv));    
    errp2m{kk} = zeros(1,numel(pv));
    
    errI1T{kk} = zeros(1,numel(Tv));
    errI2T{kk} = zeros(1,numel(Tv));
end

pvT = zeros(numel(pv),numel(Tv));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I'm going to only work with temperatures where Ts is true
Ts = true(size(Tv));

ending = zeros(size(Tv));
thisone = 0;

pvT = zeros(numel(pv),numel(Tv));

while IsayGo
    thisone = thisone + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Recompute statistics 
               %(separately for each temperature and channel, possibly using 
                    %a rough fit to the spectra gotten as a function of p)

    %Trapezoidal rule sum (interpolates linearly between points)
    % sum( (f- + f+)/2 * dp ) 

    %For integration error, Just add up discrete second derivatives and divide by 6.
    % If the sign changed a lot we could treat the parts with opposite sign as
    % independent (and therefore add in quadrature) but this should be rare
    % here anyway so we just add the absolute values.

    %Then just add the error in the values themselves in quadrature.

    % For where to sample more we just take the sum of the std and the
    % trapezoidal error at each site and we want to sample new p-values from
    % this distribution
    toc
    tic

%    II = cell(1,knum);

%     Ier = cell(1,knum); %Ina = cell(1,knum);
%     Ier1 = Ier; Ier2 = Ier;
%     Iter = cell(numel(pv),knum);
% 
%     err1 = cell(1,knum); err2 = err1;
%     errp1 = cell(1,knum);
%     errp2 = err1; errp1m = err1; errp2m = err1; errp = err1; errpm = err1;
%     errI1 = err1; errI2 = err1; errI1m = err1; errI2m = err1;

    ps = find(pvs>0);

    for kk = 2:knum

        %add them in the same loop
        Ier{kk}(Ts,:) = 0*Itt{1,kk}(Ts,:); 
        Ier1{kk}(Ts,:) = Ier{kk}(Ts,:); 
        Ier2{kk}(Ts,:) = Ier{kk}(Ts,:);
        II{kk}(Ts,:)  = 0*Itt{1,kk}(Ts,:);
        
        err1{kk}(:,Ts) = 0;
        err2{kk}(:,Ts) = 0;
        
        errp1{kk} = zeros(1,numel(pv));    
        errp2{kk} = zeros(1,numel(pv));

        errp1m{kk} = zeros(1,numel(pv));    
        errp2m{kk} = zeros(1,numel(pv));

        %There are several perspectives on error here. Probably easiest to
        %just look at the size of a variable to see what dimension it
        %retains - other variables have either been averaged or the maximum
        %has been taken over dimensions that don't remain. The total number
        %of dimensions to consider here is   Temperature, flux sector, 

        for pind = 1:numel(ps)
            %if pv(pind) > 0
                if pind ~= 1 && pind ~= numel(ps)
                    %Integration error
                    Iter{ps(pind),kk}(Ts,:) = abs( 2*Itt{ps(pind),kk}(Ts,:) - Itt{ps(pind-1),kk}(Ts,:) - Itt{ps(pind+1),kk}(Ts,:) )*( ps(pind+1) - ps(pind-1) - 2 )/12 *2/( (ps(pind)) - (ps(pind-1)) );
                    %Integral value
                    thing = ( Itt{ps(pind-1),kk}(Ts,:) + Itt{ps(pind),kk}(Ts,:) ) * ( (ps(pind)) - (ps(pind-1)) )/2;
                    if kk>1 && any(any(thing < 0))
                        error('wtf!')
                    end
                    II{kk}(Ts,:) = II{kk}(Ts,:) + thing;
                    
                    % errors where we sum over flux sectors (integrate out p)
                    %2: integral error
                    Ier2{kk}(Ts,:) = Ier2{kk}(Ts,:) + ( Iter{ps(pind),kk}(Ts,:) + Iter{ps(pind-1),kk}(Ts,:) ) * ( (ps(pind)) - (ps(pind-1)) )/2;
                    Ier1{kk}(Ts,:) = Ier1{kk}(Ts,:) + ( Its {ps(pind),kk}(Ts,:).^2 + Its {ps(pind),kk}(Ts,:).^2 ) * ( (ps(pind)) - (ps(pind-1)) )/2;
                else
                    Iter{ps(pind),kk}(Ts,:) = 0*Itt{ps(pind),kk}(Ts,:);
                    if pind == 1
                        II{kk}(Ts,:) = II{kk}(Ts,:) + Itt{ps(pind),kk}(Ts,:)/2;
                    else
                        II{kk}(Ts,:) = II{kk}(Ts,:) + ( Itt{ps(pind-1),kk}(Ts,:) + Itt{ps(pind),kk}(Ts,:) ) * ( ps(pind) - ps(pind-1) )/2;
                        II{kk}(Ts,:) = II{kk}(Ts,:) + Itt{ps(pind),kk}(Ts,:)/2;
                        
                        %2: integral error
                        Ier2{kk}(Ts,:) = Ier2{kk}(Ts,:) + ( Iter{ps(pind),kk}(Ts,:) + Iter{ps(pind-1),kk}(Ts,:) ) * ( (ps(pind)) - (ps(pind-1)) )/2;
                        %1: statistical error
                        Ier1{kk}(Ts,:) = Ier1{kk}(Ts,:) + ( Its {ps(pind),kk}(Ts,:).^2 + Its {ps(pind),kk}(Ts,:).^2 ) * ( (ps(pind)) - (ps(pind-1)) )/2;
                    end
                end

        end
        % sum variances, divide by number of independent ones
        Ier1{kk}(Ts,:) = sqrt( Ier1{kk}(Ts,:)/numel(ps) );
        %This is the error as a function of T and E
        Ier{kk}(Ts,:) = sqrt( Ier1{kk}(Ts,:).^2 + Ier2{kk}(Ts,:).^2 ); 
       
        if kk>1
            for pind = 1:numel(ps)
                %Now integrate over energy and temperature, but keep p
                
                %mean along E, then divide by average II for each T
                %(The error only matters relative to the total integral for
                %    each T)
                err1{kk}(ps(pind),Ts) = mean( Its {ps(pind),kk}(Ts,:) ,2 ) ./ mean( II{kk}(Ts,:) ,2) ;
                err2{kk}(ps(pind),Ts) = mean( Iter {ps(pind),kk}(Ts,:) ,2 ) ./ mean( II{kk}(Ts,:) ,2) ;

                %Take the mean and max along T 
                errp1{kk}(ps(pind)) = mean( err1{kk}(ps(pind),Ts) );
                errp2{kk}(ps(pind)) = mean( err2{kk}(ps(pind),Ts) );

                errp1m{kk}(ps(pind)) = max( err1{kk}(ps(pind),Ts) );
                errp2m{kk}(ps(pind)) = max( err2{kk}(ps(pind),Ts) );

                %These represent the distribution along p of where the error is
                %coming from
                errp{kk}(ps(pind)) = sqrt( errp1{kk}(ps(pind)).^2 + errp2{kk}(ps(pind)).^2 );
                errpm{kk}(ps(pind))= max(errp1m{kk}(ps(pind)),errp2m{kk}(ps(pind)));

            end
            
            %Using errors that have integrated over E and T, integrate over p
            %(leaving only a number for each kk)
            errI1{kk} = sum( interp1(ps-1,errp1{kk}(ps),0:numplaq) );
            errI2{kk} = sum( interp1(ps-1,errp2{kk}(ps),0:numplaq) );
            errI1m{kk} = sum( interp1(ps-1,errp1m{kk}(ps),0:numplaq) );
            errI2m{kk} = sum( interp1(ps-1,errp2m{kk}(ps),0:numplaq) );

            %These ones leave only T as a variable
            %errI1T{kk} = zeros(1,numel(Tv));
            %errI2T{kk} = errI1T{kk};
            for in = 1:numel(Tv)
                if Ts(in)
                    errI1T{kk}(in) = sum( interp1(ps-1,err1{kk}(ps,in),0:numplaq) );
                    errI2T{kk}(in) = sum( interp1(ps-1,err2{kk}(ps,in),0:numplaq) );
                end
            end
        end
    end

    % = mean(norma(bfs(:,ps),2),1);
    limsE = Ev;

    toc
     
    %%%%%%%%%%%%%%%%%Do error distributions for each T using err1 and
    %%%%%%%%%%%%%%%%%err2!!!
    %It's fine to force every T in a new flux sector.
    
    %% plan next step 
     
    %The error distribution along p
%    Ivar  = errpm{2}.^2+errpm{3}.^2+errpm{4}.^2+errpm{5}.^2; %errpm{1}.^2+
    Ivar  = errpm{2} + errpm{3} + errpm{4} + errpm{5}; %errpm{1}.^2+
    pvst = 2*ones(size(pvs)).'; 
    pvst(pvs ~= 0) = pvs(pvs ~= 0);
    Ivar = Ivar./pvst;
    %The distribution of expected error improvement goes like the current
    %error (std) divided by the number of times this point has been probed.
    %In case we haven't probed it yet, then the error will likely be cut in
    %half by introducing the new flux sector (in the limit that it's
    %statistical error is ignorable, otherwise just probe that sector and
    %think about its statistical error once you have a std for it).
    dist2 = interp1(ps-1,Ivar(ps),0:numplaq);

    %The error distribution along p
    subs = gendist(dist2,1,pnum);
    pv = rp*accumarray(subs.',1,[numplaq+1,1]);
       
    
    %a rough estimate of the log of the partition function for each temp.
    %Don't use a different one than your started with (pf1) or the
    %normalization of the integral will change...
    pf(Ts) = 0;
    for pind = 2:numel(ps)
        pf(Ts) = pf(Ts) + ( bfs(Ts,ps(pind)) + bfs(Ts,ps(pind-1)) )*(ps(pind)-ps(pind-1))/2;
    end
    pf(Ts) = pf(Ts) + ( bfs(Ts,ps(1)) + bfs(Ts,ps(numel(ps))) )/2;
    pf(Ts) = log(pf(Ts));
   
    
    %% compute some more stuff
    
    nup = sum(ps);
    pvals = zeros(1,nup);
    Evalys = pvals;
    for ind = 1:numel(ps)
        pind = ps(ind);
        for En = Endist{pind}
            pvals(pind) = pind-1;
            Evalys(pind)= En ;
        end
    end
    
    % mean relative error on plot
%     maxE = cell(1,5);
%     for kk = 1:5
%         maxE{kk} = max(max( Ier{kk} ))./ mean(mean( II{kk} ));
%     end
%     maxErr = max(maxE); 
%     disp(maxErr)
%     maxErrs = [maxErrs,maxErr];
%     
%     figure; plot(maxErrs)
%     xlabel('cycles');
%     ylabel('max Error');

    meanErr = zeros(1,5);
    maxMeanErr = meanErr;
    prevErr = meanErr;
    
    %meanErrFav = meanErr;
    %prevErrFav = meanErr;
    
    meanErrT = zeros(numel(Tv),5);
    
    fav = 4;
       
    if thisone == 1
        IIprev = II*0;
    end
    
    for kk = 2:5
        meanErr(kk) = sqrt(errI1{kk}.^2 + errI2{kk}.^2);
        maxMeanErr(kk) = sqrt(max(errI1T{kk}(Ts)).^2 + max(errI2T{kk}(Ts)).^2); %mean along energy, max along T
        meanErrT(:,kk) = sqrt(errI1T{kk}.^2 + errI2T{kk}.^2);
       % if kk > 1
            prevErr(kk) = max( mean(abs( norma(II{kk},2) - norma(IIprev{kk},2) ),2) ./ mean(abs(norma(II{kk},2)),2) );
            meanVal(kk) = mean(mean( norma(II{kk},2) )); %This is the mean of the plot we're after
      %  end
       % meanErrFav(kk) = meanErrT(fav,kk);
        %tickle = mean(abs( norma(II{kk},2) - norma(IIprev{kk},2) ),2) ./ mean(abs(norma(II{kk},2)),2);
        %prevErrFav(kk) = tickle(fav);
    end
    % meanErr = max(meanE);
    if thisone == 1
        meanErrs = meanErr;
        maxMeanErrs = maxMeanErr;
        prevErrs = prevErr;
        meanErrTs = {meanErrT};
        
     %   meanErrFavs = meanErrFav;
     %   prevErrFavs = prevErrFav;
    else
        meanErrs   = [meanErrs;meanErr];
        maxMeanErrs= [maxMeanErrs;maxMeanErr];
        prevErrs   = [prevErrs;prevErr];
        meanErrTs  = [meanErrTs;meanErrT];
        
     %   meanErrFavs = [meanErrFavs;meanErrFav];
      %  prevErrFavs = [prevErrFavs;prevErrFav];
    end
    disp('maxMeanErr: ')
    disp(maxMeanErr)
    disp('MeanErr: ')
    disp(meanErr)
    
    
    %Is the error low enough (all three of previous runs should be low
    %enough)
    %correct for that we are just trying to make a plot look nice, 
              %whose maxValue is set to maxVal
    
    for in = 1:numel(Tv)          
       if (thisone >= 3) 
           tem = max( [abs(meanErrTs{thisone-2}(in,2:5)) ./ meanVal(2:5), ...
               abs(meanErrTs{thisone-1}(in,2:5)) ./ meanVal(2:5),...
               abs(meanErrTs{thisone}(in,2:5)) ./ meanVal(2:5) ] );
           if tem < lowenough
               Ts(in) = false;
               ending(in) = thisone;
           end
       elseif max( abs(meanErrTs{thisone}(in,2:5)) ./ meanVal(2:5) ) < lowenough/10
           Ts(in) = false;
           ending(in) = thisone;
       end
    end
    
    disp(['Ts: ',sprintf('%d',Ts)])
     
    if ~any(Ts)
        IsayGo = false;
        break
%    else
%        Ts = true(size(Ts));
    end
    
    temp = bfs(:,ps);
    mn = mean(temp,2);
    mn(mn==0) = 1;
    bfy = bsxfun(@rdivide,temp,mn);
    
        
    %Bump up the value to at least rp everywhere that we compute to make the
    %error analysis easier/more reliable
%    new = (pv>0) & (pvs == 0); 
%     if any(new)
%         pv(new) = rp;
%         Ts = true(size(Tv));
%     end
   % disp(Ts)
    %the minimum for any sector is rp so that we always have some stats.
    pv(0<pv & pv<rp) = rp;
    
    
    %Save data after each megastep

    save(['stretch_flux2_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)]...
        ,'meanErrs','Tc','Tv','Tvis','Tcomp','Ev','II','Endist','pvs','ps','Ier','Ier1','Ier2','errp','errp1','errp2',...
        'errp1m','errp2m','errpm','bfs','bfy','prevErrs','maxMeanErrs','errI1T','errI2T','errI1',...
        'errI2','errI1m','errI2m','pf','pf1')


    ps = find(pvs~=0);
    
    %pvtot = pvtot + pv;
    %Ina,Ier,errp,errp1,errp2, Ens,Enserr,bfy ,bfs
    %EnsSig, EnSigerr, bfsErr    ,Ensdist
    
    %% Plot after each MegaStep
if plotting    
    savePlots = (mod(size(meanErrs,1),10) == 0);
%     
%     %I versus T and E (II)
%     hh=figure(h1);%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis),norma(II{5}(Tvis,:),2));
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
% 
%     % err versus T and E
%     hh=figure(h2);%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis), bsxfun( @rdivide, Ier{5}(Tvis,:), mean(II{5}(Tvis,:) ,2) ) );
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_err_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     
%         % err versus T and E
%     hh=figure(h3);%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis), bsxfun( @rdivide, Ier1{5}(Tvis,:), mean(II{5}(Tvis,:) ,2) ) );
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
% %    filename = ['stretch_err1_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
% %    saveas(hh,filename)
% %    print(hh, '-dpng', filename);
% %    print(hh, '-depsc', filename);
%     
%         % err versus T and E
%     hh=figure(h4);%('Position',position);
%     hold on;
%     uimagesc(limsE,Tv(Tvis), bsxfun( @rdivide, Ier2{5}(Tvis,:), mean( II{5}(Tvis,:),2) ) );
%     axis([0, emax, -inf, inf])
%     colorbar
%     caxis auto
%     hold off;
% %    filename = ['stretch_err2_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
% %    saveas(hh,filename)
% %    print(hh, '-dpng', filename);
%    print(hh, '-depsc', filename);

    % mean err vs p
    %DOS
%     hh=figure(h5);
%     hold on;
%     plot(ps-1,errp{5}(ps));
%     hold off;
%     filename = ['stretch_Ixy2_perr_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end

%     hh=figure(h6);
%     hold on;
%     plot(ps-1,errp1{5}(ps));
%     hold off;
%     filename = ['stretch_Ixy2_perr1_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
% 
%     hh=figure(h7);
%     hold on;
%     plot(ps-1,errp2{5}(ps));
%     hold off;
%     filename = ['stretch_Ixy2_perr2_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end
%     
%     hh=figure(h8);
%     hold on;
%     scatter(pvals(ps),Evalys(ps),5,'filled')
%     hold off;
%     %plot(ps,Ens(ps),'k',ps,Ens(ps)-Enserr(ps),'k--',ps,Ens(ps)+Enserr(ps),'k--');
%     filename = ['stretch_En_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end

%     hh=figure(h9);
%     plot(ps-1,bfy(ps));
%     xlabel('p');
%     ylabel('BF');
%     filename = ['stretch_bf_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_s_',num2str(round(1000*s))];
%     saveas(hh,filename)
%     print(hh, '-dpng', filename);
%     print(hh, '-depsc', filename);
%     
%     %I versus T and E (II) on a log plot
%     hh=figure(h13);%('Position',position);
%     hold on;
%     uimagesc(limsE,log(Tv),norma(II{2},2));
%     axis([0, emax, -5, 2])
%     colorbar
%     caxis auto
%     hold off;
%     filename = ['stretch_comp_Ixy2_finiteT_log_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
%     if savePlots
%         saveas(hh,filename)
%         print(hh, '-dpng', filename);
%         print(hh, '-depsc', filename);
%     end


    if thisone > 1
        figure(h9); hold on; plot(1./meanErrs(:,2:5)); hold off;

        figure(h10); hold on; plot(1./prevErrs(:,2:5)); hold off;

        figure(h11); hold on; plot(1./maxMeanErrs(:,2:5)); hold off;
        
       % figure(h12); hold on; plot(1./prevErrFavs(:,2:5)); hold off;
        
      %  figure(h13); hold on; plot(1./meanErrFavs(:,2:5)); hold off;
        
    end
    %Make sure I can see the plots
    pause(.01)
end
    
    
     % Manually distrubute the problem over parallel workers
    ps0 = find(pv~=0);
    q = numel(ps0);
    list = 1:q;
    parN = min(parNum,q);
    thisList = cell(1,parN);
    pvlist = cell(1,parN);
    for parInd = 1:parN
        pvlist{parInd} = zeros(size(pv));
        thisList{parInd} = ps0( mod(list,parN) == (parInd-1) );
        pvlist{parInd}( thisList{parInd} ) = pv( thisList{parInd} );
    end    
    
%     [valy, indy] = max(pv);
%     if valy > sum(pv)/2
%         %split up this run...?
%     end
    
    Tvtemp = Tv(Ts);
    pf1temp = pf1(Ts);
    %% The big step
    parfor parInd = 1:parN
        %try
        [~,IttList{parInd},ItsList{parInd},bfsList{parInd},bfserrList{parInd},EndistList{parInd}] ...
            = stretch_flux6(rmax,b,s,pvlist{parInd},Tvtemp,pf1temp,GSen); 
        %catch
        %   disp(wtf) 
        %end
    end
    
    %Unpack the distributed results from parallel eval
    Itt0 = IttList{1};
    Its0 = ItsList{1};
    bfs0 = bfsList{1};
    bfserr0 = bfsList{1};
    Endist0 = EndistList{1};
    pv2 = pvlist{1}; %I compute this just for bug checking
    for parInd = 2:parN
        for kk=1:knum
            Itt0( thisList{parInd},kk ) = IttList{parInd}( thisList{parInd},kk );
            Its0( thisList{parInd},kk ) = ItsList{parInd}( thisList{parInd},kk );
        end
        bfs0( :,thisList{parInd} ) = bfsList{parInd}( :,thisList{parInd} );
        bfserr0( :, thisList{parInd} ) = bfserrList{parInd}( :,thisList{parInd} );
        Endist0( thisList{parInd} ) = EndistList{parInd}( thisList{parInd} );
        pv2( thisList{parInd} ) = pvlist{parInd}( thisList{parInd} );
    end  
    if ~all(pv2 == pv)
        error('screwed up parallel kernel distribution')
    end
    

    tot = sum(sum(Ts));
    
    %Close the figures
    %close all;
if plotting
%     figure(h12); %14
%     plot(pv)
%     %Make sure I can see the plots
%     pause(.01)
end
    
    
    IIprev = II;
    
    
        %% incorporate this data
        
    %There is implicit use of homemade cell overides here
    %The functions are all element-wise in cell space

%There are games to play here since Inf/Inf = Nan needs to be set to 1
    %namely in the case where the other weight is 0.
    %weight = cell(1,5);
    %weight0 = cell(1,5);
    for pind = 1:numel(pv)
        if pvs(pind) > 0 && pv(pind) > 0 % thisone > 1
            % kk == 0 case
            Itt{pind,1} = ( pvs(pind)*Itt{pind,1} + pv(pind)*Itt0{pind,1} ) ...
                    / ( pvs(pind) + pv(pind) ) ;
            Its{pind,1} = sqrt( pvs(pind)*Its{pind,1}.^2 + pv(pind)*Its0{pind,1}.^2 ) ...
                    / ( pvs(pind) + pv(pind) ) ;
            %the rest
            for kk=2:5
                Itt{pind,kk}(Ts,:) = ( pvs(pind)*Itt{pind,kk}(Ts,:) + pv(pind)*reshape( Itt0{pind,kk},sum(Ts),numel(Itt0{pind,kk})/sum(Ts) ) ) ...
                    / ( pvs(pind) + pv(pind) ) ;
                Its{pind,kk}(Ts,:) = sqrt( pvs(pind)*Its{pind,kk}(Ts,:).^2 + pv(pind)*reshape( Its0{pind,kk}.^2,sum(Ts),numel(Its0{pind,kk})/sum(Ts) ) ) ...
                    / ( pvs(pind) + pv(pind) ) ;
            end
            bfs(Ts,pind) = ( pvs(pind)*bfs(Ts,pind) + pv(pind)*bfs0(:,pind) ) ...
                    / ( pvs(pind) + pv(pind) ) ;
            bfserr(Ts,pind) = sqrt( pvs(pind)*bfserr(Ts,pind).^2 + pv(pind)*bfserr0(:,pind).^2 ) ...
                / ( pvs(pind) + pv(pind) ) ;
        %end 
         elseif pv(pind)>0
            for kk = 1:5
                if numel(Itt0{pind,kk}) == 0 && numel(Its0{pind,kk}) == 0
                    warning('this is zero')
                end
                Itt{pind,kk} = zeros(numel(Tv),bins);
                Its{pind,kk} = zeros(numel(Tv),bins);
                
                Itt{pind,kk}(Ts,:) = Itt0{pind,kk} ;
                %weighted average of variances divided by sqrt(n_tot)
                Its{pind,kk}(Ts,:) = Its0{pind,kk} ;
                %end
            end            
            bfs(Ts,pind) = bfs0(:,pind) ;
            bfserr(Ts,pind) = bfserr0(:,pind);            
        end
        
        Endist{pind} = [Endist{pind},Endist0{pind}];
    end
    %Endist = cellfun(@horzcat,Endist,Endist0,'UniformOutput', false);
       
	pvs = pvs + pv;
          

end

%% Some final plots
savePlots = true;

  %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    uimagesc(limsE,Tv(Tvis),norma(II{2}(Tvis,:),2));
    axis([0, emax, -inf, inf])
    colorbar
    caxis auto
    hold off;
    filename = ['stretch_comp_DOS_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
      %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    uimagesc(limsE,Tv(Tvis),norma(II{3}(Tvis,:),2));
    axis([0, emax, -inf, inf])
    colorbar
    caxis auto
    hold off;
    filename = ['stretch_comp_Ixx_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
      %I versus T and E (II)
    hh=figure;%('Position',position);
    hold on;
    uimagesc(limsE,Tv(Tvis),norma(II{4}(Tvis,:),2));
    axis([0, emax, -inf, inf])
    colorbar
    caxis auto
    hold off;
    filename = ['stretch_comp_Ixy_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    
      %I versus T and E (II)
    hh=figure(h1);%('Position',position);
    hold on;
    uimagesc(limsE,Tv(Tvis),norma(II{5}(Tvis,:),2));
    axis([0, emax, -inf, inf])
    colorbar
    caxis auto
    hold off;
    filename = ['stretch_comp_Ixy2_finiteT_b_rmax_',num2str(round(rmax)),'_b_',num2str(round(b)),'_s_',num2str(1000*s)];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    %Energy versus p
    hh=figure(h8);
    hold on;
    scatter(pvals(ps),Evalys(ps),5,'filled')
    hold off;
    filename = ['stretch_En_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_s_',num2str(round(1000*s))];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end
    % distribution of ps samples
    hh=figure;
    hold on;
    plot(pvs)
    hold off;
    filename = ['stretch_pvs_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_s_',num2str(round(1000*s))];
    if savePlots
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
    end

stopTime = clock;
disp('Time in hours:')
disp(etime(stopTime,startTime)/3600);

diary off
