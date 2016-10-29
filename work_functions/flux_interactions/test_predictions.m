rmax = 15;
s = 0.04;
b = 10;
flag = -1; %(don't plot)
nmax = 5;

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes

load(['flux_interaction_data_n_',num2str(5),'_rmax_',num2str(rmax),'_b_',num2str(round(b)),'_s_',num2str(1000*s)])

%Get the base Hamiltonian
[H,Rxx0,Rxy0,RxyA0,RxyB0,~] = stretch_2D_makeH(rmax,b,s,0,0,flag);
En00 = -sum(svd(H));
En = 0;
%H = flip_all(rmax,H,Rxx0,Rxy0,RxyA0,RxyB0);
%En = -sum(svd(H))-En00;

%initialization
maxrels = 400;
pbool = false(numplaq,1);
pbools= false(numplaq,maxrels);
ps    = zeros(1,maxrels);
rs    = zeros(1,maxrels);
ks    = zeros(1,maxrels);
ms    = zeros(1,maxrels);
Enerrs= zeros(1,maxrels);
Enerrs0= zeros(1,maxrels);
Ens   = zeros(1,maxrels);
Ensp  = zeros(1,maxrels);

profile on

%main loop
for ii = 1:maxrels
    
    %Choose a random plaquette.
    plaq = randi(numplaq);
    eps = 0.0001; %A small number to help avoid numerical errors with ceil
    ru = ceil( ( -1+sqrt(1+(plaq-1-eps)*4/3) )/2 );
    iu = plaq - 1 - 3*ru*(ru-1);
    ku = floor( (iu - 1 + eps)/ru );
    mu = iu - ru*ku;
    if ru==0
        mu = 0;
        ku = 0;
    end
    %kk = mod(ku,2);
    ju = ru*(ru-1)/2 + mu + 1;
    
    %predict the energy cost to flipping it
    [dEnp,dEnp0] = predict_En(ru,ku,mu,pbool,Ens1,Ens2,numplaq,nmax);
    Enp = En + dEnp;
    
    %flip it
    H = flip_plaquette(ru,ku,mu,rmax,H,Rxx0,Rxy0,RxyA0,RxyB0);
    pbool(plaq) = ~pbool(plaq);
    
    %measure the energy cost to flipping it
    En1 = -sum(svd(H))-En00;
    dEnm = En1-En;
    
    %Store stuff
    pbools(:,ii) = pbool;
    ps(ii) = sum(pbool);
    rs(ii) = ru;
    ks(ii) = ku;
    ms(ii) = mu;
    
    Enerrs(ii)=dEnp-dEnm;
    Enerrs0(ii)=dEnp0-dEnm;
    Ens(ii) = dEnm;
    Ensp(ii)= dEnp;    
    
%     disp(Enerrs(ii))
%     disp(dEnp0-dEnm)
%     ps = get_NNNNs(ru,ku,mu,nmax);
%     disp(find(pbool(ps(ps<=numplaq)).'))
%     %disp(dEnm)
%     disp( Ens2(ku+1,ju, find(pbool(ps(ps<=numplaq)) ) ) )
    
    %update current energy    
    En = En1;
    
end

toc;

profile off
profile viewer

%Plot some statistics

figure;
scatter(ps,Enerrs);

figure;
scatter(ps,Enerrs0);

figure;
scatter(ps,Ens);

%figure;
%plot(1:maxrels,Enerrs);


disp('mean energy')
disp(mean(Ens))

disp('mean err')
disp(mean(Enerrs))
%disp(mean(Enerrs)/mean(Ens))

disp('std err')
%disp(std(Enerrs))
%disp(std(Enerrs)/mean(Ens))


disp(std(Enerrs))
disp(std(Enerrs0))
disp(std(Ens))

Enstds = zeros(10,1);
Enstds0 = zeros(10,1);
Enstds1 = zeros(10,1);
for v = 1:10
    Enstds(v) = std(Enerrs( (v-1)*40 + (1:40) ) );
    Enstds0(v) = std(Enerrs0( (v-1)*40 + (1:40) ) );
    Enstds1(v) = std(Ens( (v-1)*40 + (1:40) ) );
end
disp(Enstds.')
disp(Enstds0.')
disp(Enstds1.')




