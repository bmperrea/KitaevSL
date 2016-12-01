function [En,H,Es] = stretch_2D_6_en(rmax,b,s,p,~,~)
% if nargin == 5 || flag == -1
%     slow = 0;
%     if nargin == 5
%         flag = 0;
%     end
% else
%     slow = flag; %should be boolean
% end
% step = 0.03;
% 0<= rmax <~ 30 (integer number of layers) to run in a matter of seconds
% 0<= rmax <~ 50 to run in hours
% 0 < b <~ 1 (linear change in J due to change in d)
% 0 < c < 1/b is percentage strain at max and min parts of the honeycomb 
%(+y, -y , e.g.)
% p <= ( numplaq = 1 + 3*rmax*(rmax+1) ) is the number of random fluxes.

%   Here's the labeling of sites:
%             -     -
%          -     6     -
%          7     0     -
%             2     1     -
%             3     4     -
%          8     5     -

%%%%%%%(linear change in r as you go outward from center)
c = s/(sqrt(3)*rmax+sqrt(3)/2);

%Bfield = -8*c*b;
%E0 = 3*sqrt(abs(Bfield));
%disp([Bfield,E0])

%initial data
imax = 6*(rmax+1)^2;
S = imax/2;
%bins = 200;
%emax = 12.5;
%Z = imax*(emax/2)/bins;
%Ev = (1:bins)'*emax/bins;

%initialization
H = zeros(S,S);
%Rxx = sparse(S,S); % H;
%Rxy = Rxx; 
%RxyA = Rxx;
%RxyB = Rxx;

%Close the figures
% close all;
% fav_fig = figure;
% axis equal
% lw=4; lw2 =2;
% fs=26;

% DD = zeros(size(Ev)); Dd = DD;
% Ixx = DD; Iyy = DD; Ixy = DD; Ixyxx = DD; Ixyyy = DD; Ixxyy = DD;

tic

%% Build the Hamiltonian and Raman matrices

%put in the NNN bonds on the center plaquette
r0 = rr(0,0,0,0);     %even sublattice
r1 = rr(0,0,0,1);     %odd
r2 = rr(0,0,-1,1);   %odd
r3 = rr(-1,0,-1,0); %even
% r4 = rr(-1,0,0,0);   %even
% r5 = rr(-1,0,-1,1);   %odd

i0 = 3;
i1 = 1;
i2 = 3;
i3 = 2; 
% i4 = 1;
% i5 = 2;

% RxyA(i0,i4) = 1; draw_line2(r0,r4)
% RxyA(i4,i3) = 1; draw_line2(r4,r3)
% RxyA(i3,i0) = 1; draw_line2(r3,r0)
% RxyB(i1,i5) = 1; draw_line2(r1,r5)
% RxyB(i5,i2) = 1; draw_line2(r5,r2)
% RxyB(i2,i1) = 1; draw_line2(r2,r1)

for r=0:rmax
    for k=0:5
        
        r1 = rr(r,k,0,1);
        l  = (2*r+1)*k;
        i1 = ind(r,l+1);
        
        %Put three bonds on each plaquette
        for m=1:r
    l = (2*r+1)*k + 2*m;
                          %for k even this is on the
    r0 = rr(r,k,m,0);     %even sublattice
    r1 = rr(r,k,m,1);     %odd
    r2 = rr(r,k,m-1,1);   %odd
    r3 = rr(r-1,k,m-1,0); %even
%    r4 = rr(r-1,k,m,0);   %even
%    r5 = rr(r-1,k,m-1,1);   %odd
    
    i0 = ind(r,l);
    i1 = ind(r,l+1);
    i2 = ind(r,l-1);
    i3 = ind(r-1,l-2-2*k); 
%    if m<r
 %       i4 = ind(r-1,l-2*k);
%    else
  %      i4 = ind(r,l+2);
%    end
 %   i5 = ind(r-1,l-1-2*k);
    
%    draw_text(r5,i5)
    
%     if r==rmax
%         disp( dd(r0,r1,c) );
%     end
    
    if mod(k,2)==0
        H  (i0,i1) = jj(r0,r1);
%         Rxx(i0,i1) = jj(r0,r1)*dd(r0,r1,1)^2;
%         Rxy(i0,i1) = jj(r0,r1)*dd(r0,r1,1)*dd(r0,r1,2); 
%         draw_line(r0,r1)
%         draw_text(r0,i0)
%         draw_text(r1,i1)
        
        H  (i0,i2)   = jj(r0,r2);
%         Rxx(i0,i2) = jj(r0,r2)*dd(r0,r2,1)^2;       
%         Rxy(i0,i2) = jj(r0,r2)*dd(r0,r2,1)*dd(r0,r2,2);
%         draw_line(r0,r2)
%         draw_text(r2,i2)
%         draw_text(r4,i4)
        
        if k == 0 && m == 1
            lmax = 6*(r^2 - (r-1)^2);
            i3 = ind(r-1,lmax);
%         elseif k == 5 && m == r
%             lmin = 2;
%             i4 = ind(r-1,lmin);
        end

        H  (i3,i2) = jj(r3,r2);
%         Rxx(i3,i2) = jj(r3,r2)*dd(r3,r2,1)^2;
%         Rxy(i3,i2) = jj(r3,r2)*dd(r3,r2,1)*dd(r3,r2,2);
%         draw_line(r3,r2)
%         draw_text(r3,i3)
        
%         RxyA(i0,i4) = 1; draw_line2(r0,r4)
%         RxyA(i4,i3) = 1; draw_line2(r4,r3)
%         RxyA(i3,i0) = 1; draw_line2(r3,r0)
%         RxyB(i1,i5) = 1; draw_line2(r1,r5)
%         RxyB(i5,i2) = 1; draw_line2(r5,r2)
%         RxyB(i2,i1) = 1; draw_line2(r2,r1)
        
    else %k is odd        
        
 %       if k==5 && m==r %last plaquette
 %           lmin = 2;
 %           i4 = ind(r,lmin-1);
 %       end
%        draw_text(r4,i4)
        
        H  (i1,i0) = jj(r0,r1);
%         Rxx(i1,i0) = jj(r0,r1)*dd(r0,r1,1)^2;
%         Rxy(i1,i0) = jj(r0,r1)*dd(r0,r1,1)*dd(r0,r1,2);
%         draw_line(r1,r0)
%         draw_text(r0,i0)
%         draw_text(r1,i1)
        
        H  (i2,i0) = jj(r0,r2);
%         Rxx(i2,i0) = jj(r0,r2)*dd(r0,r2,1)^2;
%         Rxy(i2,i0) = jj(r0,r2)*dd(r0,r2,1)*dd(r0,r2,2);
%         draw_line(r2,r0)
%         draw_text(r2,i2)
        
        H  (i2,i3) = jj(r3,r2);
%         Rxx(i2,i3) = jj(r3,r2)*dd(r3,r2,1)^2;
%         Rxy(i2,i3) = jj(r3,r2)*dd(r3,r2,1)*dd(r3,r2,2); 
%         draw_line(r2,r3)
%         draw_text(r3,i3)
        
%         RxyB(i0,i4) = 1; draw_line2(r0,r4)
%         RxyB(i4,i3) = 1; draw_line2(r4,r3)
%         RxyB(i3,i0) = 1; draw_line2(r3,r0)
%         RxyA(i1,i5) = 1; draw_line2(r1,r5)
%         RxyA(i5,i2) = 1; draw_line2(r5,r2)
%         RxyA(i2,i1) = 1; draw_line2(r2,r1)
    end
     
        end
               
        %Put one bond on every corner
        corner = { [0,-1] ,[-sqrt(3),-1]/2 , [-sqrt(3),1]/2 , [0,1] , [sqrt(3),1]/2 , [sqrt(3),-1]/2}; 
        r4 = r1 + corner{k+1};
          
        i4 = ind(r,l+2);
        %if k == 4 && r==0
        %   i4 = ind(r,1);
        if k==5
           %i1 = ind(r,0);
           i4 = ind(r,1);
        end

%        disp([i1,i4])
        
        if mod(k,2)==0
            H(  i4,i1) = jj(r4,r1);   
%             Rxx(i4,i1) = jj(r4,r1)*dd(r4,r1,1)^2; 
%             Rxy(i4,i1) = jj(r4,r1)*dd(r4,r1,1)*dd(r4,r1,2);
%             draw_line(r1,r4)
%             draw_text(r4,i4)
%             draw_text(r1,i1)
        else %k is odd
            H(  i1,i4) = jj(r4,r1);   
%             Rxx(i1,i4) = jj(r4,r1)*dd(r4,r1,1)^2; 
%             Rxy(i1,i4) = jj(r4,r1)*dd(r4,r1,1)*dd(r4,r1,2);
%             draw_line(r4,r1)
%             draw_text(r4,i4)
%             draw_text(r1,i1)
        end
    end 
end

%% Put in fluxes

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes
%If there's more fluxes asked for than plaquettes reduce to max
if p > numplaq
    disp(['p too large, using max fluxes: ',num2str(numplaq)])
    p = numplaq;
end
%If we are setting more than half then change the base gauge
if p > numplaq/2
    %set the base gauge to all fluxes
    for r=1:rmax
        for k=0:5
            for m = 1:r
                %Change the (2,3) bond and the associated NNN bonds
                
                l = (2*r+1)*k + 2*m;
%     i0 = ind(r,l);
%     i1 = ind(r,l+1);
%     i2 = ind(r,l-1);
%     i3 = ind(r-1,l-2-2*k); 
%     i4 = ind(r-1,l-2*k);
%     i5 = ind(r-1,l-1-2*k);
                i0 = ind(r,l);
               % i1 = ind(r,l+1);
                i2 = ind(r,l-1); 
                i3 = ind(r-1,l-2-2*k);
                %i5 = ind(r-1,l-1-2*k);
              %  i7 = ind(r,l-2);
               % i8 = ind(r-1,l-3-2*k);
                if k == 0 && m == 1
                    lmax = 6*(r^2 - (r-1)^2);
                    i3 = ind(r-1,lmax);
              %      i8 = ind(r-1,lmax-1);
                    lmax = 6*((r+1)^2 - r^2);
                    %i7 = ind(r,lmax);
                end
            %    if [r,k,m] == [1,1,1] %a random subcase
             %       lmax = 6*(r^2 - (r-1)^2);
               %     i8 = ind(r-1,lmax);
             %   end
              
                %Just need positions for drawing fluxes in the lattice pic
                r0 = rr(r,k,m,0);    % draw_text(r0,i0);
 %               r1 = rr(r,k,m,1);     draw_text(r1,i1);
                r2 = rr(r,k,m-1,1);  % draw_text(r2,i2);
                r3 = rr(r-1,k,m-1,0);% draw_text(r3,i3);
               % r4 = rr(r-1,k,m,0);  draw_text(r0,i0);
                %r5 = rr(r-1,k,m-1,1); %  draw_text(r5,i5);
                %r7 = rr(r,k,m-1,0);  % draw_text(r7,i7);
              %  r8 = rr(r-1,k,m-2,1); %draw_text(r8,i8);              

                %Change the bond for every other plaquette 
                    %going around each ring level
                sn  = -(-1)^m;
                sn2 = -(-1)^r * sn;    %switch for odd/even k, but only if 
                                        % r is odd.
                if mod(k,2) == 0
                    H  (i3,i2) = sn*H  (i3,i2);
%                     Rxx(i3,i2) = sn*Rxx(i3,i2);
%                     Rxy(i3,i2) = sn*Rxy(i3,i2);
%                     draw_line(r3,r2,H(i3,i2))
                  
%                     RxyA(i3,i0) = sn*RxyA(i3,i0);
%                     draw_line2(r3,r0,RxyA(i3,i0))
%                     RxyB(i5,i2) = sn*RxyB(i5,i2);
%                     draw_line2(r5,r2,RxyB(i5,i2))
%                     RxyA(i7,i3) = sn*RxyA(i7,i3);
%                     draw_line2(r7,r3,RxyA(i7,i3))
%                     RxyB(i2,i8) = sn*RxyB(i2,i8);
%                     draw_line2(r2,r8,RxyB(i2,i8))
                else
                    H  (i2,i3) = -sn2*H  (i2,i3);
%                     Rxx(i2,i3) = -sn2*Rxx(i2,i3);
%                     Rxy(i2,i3) = -sn2*Rxy(i2,i3);
%                     draw_line(r2,r3,H(i2,i3))
% 
%                     RxyB(i3,i0) = -sn2*RxyB(i3,i0);
%                     draw_line2(r3,r0,RxyB(i3,i0))
%                     RxyA(i5,i2) = -sn2*RxyA(i5,i2);
%                     draw_line2(r5,r2,RxyA(i5,i2))
%                     RxyB(i7,i3) = -sn2*RxyB(i7,i3);
%                     draw_line2(r7,r3,RxyB(i7,i3))
%                     RxyA(i2,i8) = -sn2*RxyA(i2,i8); 
%                     draw_line2(r2,r8,RxyA(i2,i8))                       
                end
                              
            end
        end
    end
    %set the center plaquette with the gauge string ending at the k=2 site
    %This amounts to changing the z-bond on the right side of the plaq.
    add_flux(0,0,1)
    
    %set p to be the number of random non-flux sites
    p = numplaq - p;
end

%Generate a random flux pattern with p fluxes
%ms = []; ks = ms; rs = ms;
inds = [];
plaqs = 3*rmax*(rmax+1)+1;
while numel(inds) < p
    rans = randi(plaqs,p-numel(inds),1);
    inds = [inds;rans];
    %Remove duplicates
    inds = unique(inds);
    
%     for pp = 1:(p-numel(inds))
%         r = randi(rmax+1)-1;
%         k = randi(6)-1;
%         if r>0
%             m = randi(r);
%         else
%             m = 0;
%         end
%         rs = [rs,r];  
%         ms = [ms,m];
%         ks = [ks,k];
%         inds = [inds,ind(r,k*(2*r+1)+l)];
%     end
    %inds = C; % = inds(ia)
    %ks = ks(ia);
   % ms = ms(ia);
  %  disp( numel(inds) )
end
% if numel(ls)<p
%     error(['ls_',num2str(numel(ls))),'_p_',num2tr(p)])
% end

%Put in fluxes
if p ~= numel(inds);
    disp('the number of fluxes is wrong')
end
eps = 0.0001; %A small number to help avoid numerical errors with ceil
for pp= 1:p
    %recover (r,k,m) from the index of the plaquette
    plaq = inds(pp);
    ru = ceil( ( -1+sqrt(1+(plaq-1-eps)*4/3) )/2 );
    iu = plaq - 1 - 3*ru*(ru-1);
    ku = floor( (iu - 1 + eps)/ru );
    mu = iu - ru*ku;
    
    %Gauge this with a string of gauge transformations coming from the
    %edge found by shooting a ray in the k-direction
    add_flux(ru,ku,mu)
%     for r1 = r:rmax
%         l = 2*m + k*(2*r1+1);
%         
%         %grab the site indices on the bond to change
%         i0 = ind(r1,l);
%         i1 = ind(r,l+1);
%         i2 = ind(r,l-1);
%         i3 = ind(r-1,l-2-2*k); 
%         
%         %Change the bond
%         if mod(k,2)==0 %Then the 0 site is an even site                
%             H(i2,i0) = -H(i2,i0);
%             Rxx(i2,i0) = -Rxx(i2,i0);
%             Rxy(i2,i0) = -Rxy(i2,i0);
%             draw_line(r2,r0,H(i2,i0))
%             
%             if k == 0 && m == 1
%                 lmax = 6*(r^2 - (r-1)^2);
%                 i3 = ind(r-1,lmax);
%             end
%        
%             RxyA(i3,i0) = -RxyA(i3,i0);
%             RxyB(i2,i1) = -RxyB(i2,i1);
%         else %0 site is odd
%             H(i0,i2) = -H(i0,i2);
%             Rxx(i0,i2) = -Rxx(i0,i2);
%             Rxy(i0,i2) = -Rxy(i0,i2); 
%             draw_line(r3,r2,H(i0,i2))
%        
%             RxyB(i3,i0) = -RxyB(i3,i0);
%             RxyA(i2,i1) = -RxyA(i2,i1);          
%         end
        % update i0 for the next row
       % i0 = i0 + 2+2*k;
 %   end
end

%add_flux(1,0,1)

% RxyA = RxyA - RxyA';
% RxyB = RxyB - RxyB';

%% Compute and plot spectra

%toc
%tic

%memory

%Rxx = sparse(Rxx);
%Ryy = sparse(Ryy);
%Rxy = sparse(Rxy);

toc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonalizing
% % H2 = full(H);
% 
% %[U,D,V] = svd(H);
% memory
% clearvars H
% 
% toc
% tic
% 
% en2 = diag(D);
% clearvars D
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Histogram it
% %Rxx = U'*Rxx*U;
% Rxx = U'*Rxx*V;
% Ryy = U'*Ryy*V;
% Rxy = U'*Rxy*V;
% 
% memory
% clearvars U
% clearvars V
% 
% Rxx = (Rxx' - Rxx)/2;
% Ryy = (Ryy' - Ryy)/2;
% Rxy = (Rxy' - Rxy)/2;
% 
% %Store the energies in a big matrix
% ent = repmat(en2,1,S) + repmat(en2',S,1);
% 
% ins = (emax/2>=ent)&(ent>=0);
% 
% Wxx = 2*pi*real( conj(Rxx).*Rxx );
% [histw, histv] = histwv(2*ent(ins),Wxx(ins),0,emax,bins);
% DD  = histv/Z;
% Ixx = histw/Z;
% clearvars Wxx
% 
% Wyy = 2*pi*real( conj(Ryy).*Ryy );
% [histw, ~    ] = histwv(2*ent(ins),Wyy(ins),0,emax,bins);
% Iyy = Iyy + histw/Z;
% clearvars Wyy
% 
% Wxy = 2*pi*real( conj(Rxy).*Rxy );
% [histw, ~    ] = histwv(2*ent(ins),Wxy(ins),0,emax,bins);
% Ixy = Ixy + histw/Z;
% clearvars Wxy
% 
% Wxxyy = pi*real( conj(Rxx).*Ryy + conj(Ryy).*Rxx );
% [histw, ~    ] = histwv(2*ent(ins),Wxxyy(ins),0,emax,bins);
% Ixxyy = Ixxyy + histw/Z;
% clearvars Wxxyy 
% 
% Wxyxx = pi*real( conj(Rxx).*Rxy + conj(Rxy).*Rxx );
% [histw, ~    ] = histwv(2*ent(ins),Wxyxx(ins),0,emax,bins);
% Ixyxx = Ixyxx + histw/Z;
% clearvars Wxyxx
% 
% Wxyyy = pi*real( conj(Rxy).*Ryy + conj(Ryy).*Rxy );
% [histw, ~    ] = histwv(2*ent(ins),Wxyyy(ins),0,emax,bins);
% Ixyyy = Ixyyy + histw/Z;
% clearvars Wxyyy
% 
% %Store the 1p-DOS as well
% ins2 = (emax/4>=en2) & en2>=0;
% [~, histv] = histwv(2*en2(ins2),0*en2(ins2),0,emax/2,bins);
% Dd = Dd + 2*histv/Z;
% 
% %display the number of points above emax
% disp( sum( ceil( abs( en2(en2>emax/4) ) ) ) / ( size(en2,1)*size(en2,2) ) )
% disp( sum( ceil( abs( ent(ent>emax/2) ) ) ) / ( size(ent,1)*size(ent,2) ))

%try 
Es = svd(H);
En = -sum(Es);
  %  [En,Ev,Ddd,Dd,Ixx,Ixy,Ixy2] = dos1_loop(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax);
%catch
%    addpath('C:\Users\Brent\SkyDrive\Matlab\Monte_carlo')
%    [En,Ev,Ddd,Dd,Ixx,Ixy,Ixy2] = dos1(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax);
%end

toc

% I = cell(1,6);
% I{1} = Ev;    I{2} = Ddd;
% I{3} = Ixx;   I{4} = Ixy;
% I{5} = Ixy2;  I{6} = Dd;
% 
% %This flag is for if we don't want to make intermediate plots
% if flag ~= -1
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Plot stuff
% 
%     %initial data
%     width = 5.1;     % Width in inches
%     height = 3;    % Height in inches
%     %alw = 1;       % AxesLineWidth
%     fsz = 15;      % Fontsize
%     %fszd = 22;
%     fna = 'Helvetica'; %Fontname
%     lw = 1.5;      % LineWidth
%     msz = 8;       % MarkerSize
%     interp = 'tex';
% 
%     % The properties we've been using in the figures
%     set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
%     set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
%     set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
%     set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% 
%     set(0,'defaultAxesFontName',fna);
%     set(0,'defaultAxesFontSize',fsz);
%     set(0,'defaultTextInterpreter',interp);
% 
%     % Set the default Size for display
%     defpos = get(0,'defaultFigurePosition');
%     set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
% 
%     % Set the defaults for saving/printing to a file
%     set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
%     set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
%     defsize = get(gcf, 'PaperSize');
%     left = (defsize(1)- width)/2;
%     bottom = (defsize(2)- height)/2;
%     defsize = [left, bottom, width, height];
%     set(0, 'defaultFigurePaperPosition', defsize);
% 
%     %Save the honeycomb figure
%     if rmax<20
%         filename = ['stretch_fig_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%             '_c_',num2str(round(1000*s)),'_p_',num2str(p)];
%         saveas(fav_fig,filename)
%         print(fav_fig, '-dpng', filename);
%         print(fav_fig, '-depsc', filename);
%     end
% 
%     peaks = 7;
% 
%     %Plot DOS
%     hh=figure;
%     hold on;
%     plot(Ev/2,Ddd);
%     xlim([0,6.45])
%     %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%     %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
%     xlabel('\omega/J');
%     ylabel('DOS');
%     %set(gca,'XTick',-3:3); 
%     %set(gca,'YTick',2*(0:5));
%     %hold off;
%     yl = ylim;
%     for n = 0:peaks
%         %[Xi, ~] = ds2nfu(sqrt(n)*E0, 0);
%         %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')    
%         line(sqrt(n)*[E0,E0],yl,'Color',[.7,.7,.7])
%     end
%     plot(Ev/2,Ddd);
%     hold off;
%     filename = ['stretch_DOS_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_c_',num2str(round(1000*s)),'_p_',num2str(p)];
%     saveas(hh,filename)
%     print(hh, '-dpng', filename);
%     print(hh, '-depsc', filename);
% 
%     %Plot the ramans
%     hh=figure;%('Position',position);
%     hold on;
%     plot(Ev,Ixy2);
%     xlim([0,12.9])
%     %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%     %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
%     xlabel('\omega/J');
%     ylabel('I(\omega)');
%     legend({'I_{[xy]}'}, 'Location', 'NorthEast');
%     yl = ylim;
%     for n = 0:peaks
%         %[Xi, ~] = ds2nfu(2*sqrt(n)*E0, 0);
%         %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1') 
%         line(2*sqrt(n)*[E0,E0],yl,'Color',[.7,.7,.7])
%     end
%     plot(Ev,Ixy2);
%     hold off;
%     filename = ['stretch_I3_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_c_',num2str(round(1000*s)),'_p_',num2str(p)];
%     saveas(hh,filename)
%     print(hh, '-dpng', filename);
%     print(hh, '-depsc', filename);
% 
%     hh=figure;%('Position',position);
%     hold on;
%     plot(Ev,Ixx,Ev,Ixy);
%     xlim([0,12.9])
%     %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%     %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
%     xlabel('\omega/J');
%     ylabel('I(\omega)');
%     legend({'I_{xx}', 'I_{xy}'}, 'Location', 'NorthEast');
%     yl = ylim;
%     for n = 0:peaks-1
%         %[Xi, ~] = ds2nfu((sqrt(n)+sqrt(n+1))*E0, 0);
%         %annotation('textarrow',[Xi,Xi],[0.9,.8])%,'String','1')  
%         line((sqrt(n) + sqrt(n+1))*[E0,E0],yl,'Color',[.7,.7,.7])
%     end
%     plot(Ev,Ixx,Ev,Ixy);
%     hold off;
%     filename = ['stretch_I2_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_c_',num2str(round(1000*s)),'_p_',num2str(p)];
%     saveas(hh,filename)
%     print(hh, '-dpng', filename);
%     print(hh, '-depsc', filename);
% 
% 
%     hh=figure;%('Position',position);
%     hold on;
%     plot(Ev,Ixx,Ev,Ixy2/10,Ev/2,Ddd,Ev,Ddd);
%     %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%     %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(jx),', h=',num2str(h)])
%     xlabel('\omega/J');
%     ylabel('I(\omega)');
%     legend({'I_{xx}', 'I_{[xy]}/10','DOS(\omega)','DOS(2\omega)'}, 'Location', 'NorthEast');
%     hold off;
%     filename = ['stretch_Imix_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
%         '_c_',num2str(round(1000*s)),'_p_',num2str(p)];
%     saveas(hh,filename)
%     print(hh, '-dpng', filename);
%     print(hh, '-depsc', filename);
%     
%end


%% Nested functions

    function add_flux(rrc,k,m)  
        %This function edits H and R
        if rrc==0
            add_flux(1,5,1);
            
            %Still need to change the top-left bond of 0-plaquette
            i0 = 3;
            i2 = 3;
            i1 = 1;
            i3 = 2;
  %          i6 = 4;
            %i7 = 11;
            
            %These are only needed to draw the lattice picture!
            %Just need positions for drawing fluxes in the lattice pic
            r0 = rr(0,0,0,0);   %  draw_text(r0,i0)
            r1 = rr(0,0,0,1);   %  draw_text(r1,i1)
            r2 = rr(0,0,-1,1);  %  draw_text(r2,i2)
            r3 = rr(-1,0,-1,0); %  draw_text(r3,i3)
            %r4 = rr(r-1,k,m,0);  
            %r5 = rr(r-1,k,m,1);   %odd
       %     rc = 0;
      %      if rmax>rc
  %              r6 = rr(1,0,0,1);   %  draw_text(r6,i6)
        %        r7 = rr(0,0,-1,0);  %  draw_text(r7,i7)
    %        end
            %r8 = rr(r-1,k,m-1,1);  
            
            %rc = 0;
            flip_0_2_bond()  
            
        else

            for rc = rrc:rmax    % shoot a gauge string to the edge
                %Flip the (0,2) bond and the corresponding NNN bonds
                l = (2*rc+1)*k+2*m;
                i0 = ind(rc,l);
                i1 = ind(rc,l+1);
                i2 = ind(rc,l-1);
                i3 = ind(rc-1,l-2-2*k);
                if k == 0 
                    lmax = 6*(rc^2 - (rc-1)^2);
                    i3 = ind(rc-1,lmax);
                end

                %These are only needed to draw the lattice picture!
                %Just need positions for drawing fluxes in the lattice pic
    %            disp([r,k,m])
                r0 = rr(rc,k,m,0);    % draw_text(r0,i0)
                r1 = rr(rc,k,m,1);    % draw_text(r1,i1)
                r2 = rr(rc,k,m-1,1);  % draw_text(r2,i2)
                r3 = rr(rc-1,k,m-1,0);% draw_text(r3,i3)
                %r4 = rr(r-1,k,m,0);  draw_text(r4,i4)
                %r5 = rr(r-1,k,m,1);   draw_text(r5,i5)

            %    if rc < rmax
     %               i6 = ind(rc+1,l+2*k+1);
              %      i7 = ind(rc,l-2);
   %                 r6 = rr(rc+1,k,m,1); %  draw_text(r6,i6)
             %       r7 = rr(rc,k,m-1,0); %  draw_text(r7,i7)
                    %r8 = rr(r-1,k,m-1,1);   
           %     end

                flip_0_2_bond()

            end    
        end
        
        function flip_0_2_bond()

            %Change everything originating from the (2,0) bond!
            %(0,2)
            if mod(k,2) == 0
                H  (i0,i2) = -H  (i0,i2);
%                 Rxx(i0,i2) = -Rxx(i0,i2);
%                 Rxy(i0,i2) = -Rxy(i0,i2);  
%                 draw_line(r0,r2,H(i0,i2));
%                 RxyA(i3,i0) = -RxyA(i3,i0);
%                 draw_line2(r3,r0,RxyA(i3,i0))
%                 RxyB(i2,i1) = -RxyB(i2,i1); 
%                 draw_line2(r2,r1,RxyB(i2,i1))
              %  if rc < rmax
                    %disp([rc,rmax])
%                     RxyA(i0,i7) = -RxyA(i0,i7);
%                     draw_line2(r0,r7,RxyA(i0,i7))
%                     RxyB(i6,i2) = -RxyB(i6,i2); 
%                     draw_line2(r6,r2,RxyB(i6,i2))     
              %  end
            else
                H  (i2,i0) = -H  (i2,i0);
%                 Rxx(i2,i0) = -Rxx(i2,i0);
%                 Rxy(i2,i0) = -Rxy(i2,i0);  
%                 draw_line(r2,r0,H(i2,i0)); 
% 
%                 RxyB(i3,i0) = -RxyB(i3,i0);
%                 draw_line2(r3,r0,RxyB(i3,i0))
%                 RxyA(i2,i1) = -RxyA(i2,i1); 
%                 draw_line2(r2,r1,RxyA(i2,i1))
%                 if rc < rmax
%                     RxyB(i0,i7) = -RxyB(i0,i7);
%                     draw_line2(r0,r7,RxyB(i0,i7))
%                     RxyA(i6,i2) = -RxyA(i6,i2); 
%                     draw_line2(r6,r2,RxyA(i6,i2))   
%                 end
            end
        
        end
        
    end

    function in = ind(r,l)
        %if nargin == 2
            in = 6*r^2 + l;
       % end
        %This puts 1 -> 1 , 2->1 , 3->2 , 4->2 ,...
        in = floor( (in+1)/2 );
%         if in == 13
%             m = fk;
%         end
    end

    function r1 = rr(r,k,m,l1)
        interc = { [0,0]          , [sqrt(3),-1]/2, [sqrt(3),-3]/2 , [0,-2]         , [-sqrt(3),-3]/2 , [-sqrt(3),-1]/2 };
        sloper = { [-sqrt(3),3]/2 , [sqrt(3),3]/2 , [sqrt(3),0]    , [sqrt(3),-3]/2 , [-sqrt(3),-3]/2 , [-sqrt(3),0]    };
        slopem = { [sqrt(3),0]    , [sqrt(3),-3]/2, [-sqrt(3),-3]/2, [-sqrt(3),0]   , [-sqrt(3),3]/2  , [sqrt(3),3]/2   };
        slopel = { [sqrt(3),-1]/2 , [0,-1]        , [-sqrt(3),-1]/2, [-sqrt(3),1]/2 , [0,1]           , [sqrt(3),1]/2   };
        r1 = interc{k+1} + r*sloper{k+1} + m*slopem{k+1} + (l1)*slopel{k+1} + [0,1];
    end

    function u = uu(r1)
        vec = r1;
        x = vec(1); y = vec(2);
        u = [2*x*y,x^2-y^2];
    end

    function j = jj(r1,r2)
        j = 1 + b - b*norm( dd(r1,r2) );      
        if j <= 0
            error('j is gone!')
        end
    end

    function d = dd(r1,r2,in)
    %The third argument is optional
        du = ( uu(r1) - uu(r2) );
        d1 = r1 - r2 + c*du;
        if nargin == 3
            d = d1(in);
        else
            d = d1;
        end        
        if c*norm(du) > 1
            error('too much strain')
        end
    end

%     function draw_line(p1,p2,h)
%         p1 = p1 + c*uu(p1);
%         p2 = p2 + c*uu(p2);
%         if rmax < 12 && flag ~= -1
%             if slow
%                 pause(step)
%             end
%             if nargin == 2 || h > 0
%                 col = [0 0 .7]; %blue
%             else
%                     if h == 0
%                         dbstack
%                     end
%                 col = [.7 0 0]; %red
%             end
%             line([p1(1),p2(1)],[p1(2),p2(2)],'LineWidth',lw,'Color',col)
%         end
%     end
% 
%     function draw_line2(p1,p2,h)
%         p1 = p1 + c*uu(p1);
%         p2 = p2 + c*uu(p2);
%         if rmax < 5 && flag ~= -1
%             if slow
%                 pause(step)
%             end
%             if nargin == 2 || h > 0
%                 col = [0.3 0.3 .7]; %blue
%             else
%                     if h == 0
%                         dbstack
%                     end
%                 col = [.7 0.3 0.3]; %red
%             end
%             line([p1(1),p2(1)],[p1(2),p2(2)],'LineWidth',lw2,'Color',col,'LineStyle','--')
%         end
%     end
% 
%     function draw_text(p1,txt)
%         p1 = p1 + c*uu(p1);
%         if rmax < 5 && flag ~= -1
%             if slow
%                 pause(step)
%             end
%             text(p1(1),p1(2),num2str(txt),'FontSize',fs)
%         end
%     end


end