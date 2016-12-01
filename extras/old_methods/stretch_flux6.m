function [Ev,Itt,Its,bfs,bfserr,Endist,pv] ...
    = stretch_flux6(rmax,b,s,pv,Tv,pf1,GSen)
%bfy,Ier,Ier1,Ier2,errp1,errp2,errp1m,errp2m,errp,errpm,

%rmax = 20; rp =10; b=10; s=.04; pv~100 entries from 1 to numplaq
%Tv ~10 to 100 entries 0 to inf 

% 0<= rmax <~ 30 (integer number of layers) to run in a matter of seconds
% 0<= rmax <~ 50 to run in hours
% 0 < b <~ 1 (linear change in J due to change in d)
% 0 < c < 1/b is percentage strain at max and min parts of the honeycomb 
%(+y, -y , e.g.)
% p <= ( numplaq = 1 + 3*rmax*(rmax+1) ) is the number of random fluxes.
    %pv is the vector of ps
% Tv is a vector of temps

%   Here's the labeling of sites:
%             -     -
%          -     6     -
%          7     0     -
%             2     1     -
%             3     4     -
%          8     5     -

%%%%%%%(linear change in r as you go outward from center)
c = s/(sqrt(3)*rmax+sqrt(3)/2);

Bfield = -4*c*b;
E0 = 3*sqrt(2*abs(Bfield));
disp([Bfield,E0])

%initial data
imax = 6*(rmax+1)^2;
S = imax/2;
bins = 200;
emax = 12.0;
%Z = imax*(emax/2)/bins;
%Ev = (1:bins)'*emax/bins;

%initialization
H = zeros(S,S);
Rxx = sparse(S,S); % H;
Rxy = Rxx; 
RxyA = Rxx;
RxyB = Rxx;


% DD = zeros(size(Ev)); Dd = DD;
% Ixx = DD; Iyy = DD; Ixy = DD; Ixyxx = DD; Ixyyy = DD; Ixxyy = DD;
knum = 5;
rp = max(pv);
%{En,Ddd,Ixx,Ixy,Ixy2}
It  = cell(rp,knum);  

Itt  = cell(numel(pv),knum);  
Its  = cell(numel(pv),knum);   
%Iter1 = cell(numel(pv),5);
%Iter2 = cell(numel(pv),5);
%Iter = cell(numel(pv),knum);

%bfst = cell( rp, numel(pv) );
bfs = zeros( numel(Tv), numel(pv) );
bfserr = bfs;
tmp = zeros(numel(Tv) ,rp);
%Etmp = zeros(1,rp);
Endist = cell(1,numel(pv));
%Ens = zeros(1,numel(pv));
%Enserr = zeros(1,numel(pv));


tic

%% Build the Hamiltonian and Raman matrices

%put in the NNN bonds on the center plaquette
r0 = rr(0,0,0,0);     %even sublattice
r1 = rr(0,0,0,1);     %odd
r2 = rr(0,0,-1,1);   %odd
r3 = rr(-1,0,-1,0); %even
r4 = rr(-1,0,0,0);   %even
r5 = rr(-1,0,-1,1);   %odd

i0 = 3;
i1 = 1;
i2 = 3;
i3 = 2; 
i4 = 1;
i5 = 2;

RxyA(i0,i4) = 1; 
RxyA(i4,i3) = 1; 
RxyA(i3,i0) = 1; 
RxyB(i1,i5) = 1; 
RxyB(i5,i2) = 1; 
RxyB(i2,i1) = 1; 

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
    r4 = rr(r-1,k,m,0);   %even
    r5 = rr(r-1,k,m-1,1);   %odd
    
    i0 = ind(r,l);
    i1 = ind(r,l+1);
    i2 = ind(r,l-1);
    i3 = ind(r-1,l-2-2*k); 
    if m<r
        i4 = ind(r-1,l-2*k);
    else
        i4 = ind(r,l+2);
    end
    i5 = ind(r-1,l-1-2*k);
    
    if mod(k,2)==0
        H  (i0,i1) = jj(r0,r1);
        Rxx(i0,i1) = jj(r0,r1)*dd(r0,r1,1)^2;
        Rxy(i0,i1) = jj(r0,r1)*dd(r0,r1,1)*dd(r0,r1,2); 
        
        H  (i0,i2)   = jj(r0,r2);
        Rxx(i0,i2) = jj(r0,r2)*dd(r0,r2,1)^2;       
        Rxy(i0,i2) = jj(r0,r2)*dd(r0,r2,1)*dd(r0,r2,2);
        
        if k == 0 && m == 1
            lmax = 6*(r^2 - (r-1)^2);
            i3 = ind(r-1,lmax);
%         elseif k == 5 && m == r
%             lmin = 2;
%             i4 = ind(r-1,lmin);
        end

        H  (i3,i2) = jj(r3,r2);
        Rxx(i3,i2) = jj(r3,r2)*dd(r3,r2,1)^2;
        Rxy(i3,i2) = jj(r3,r2)*dd(r3,r2,1)*dd(r3,r2,2);
        
        RxyA(i0,i4) = 1; 
        RxyA(i4,i3) = 1; 
        RxyA(i3,i0) = 1; 
        RxyB(i1,i5) = 1;
        RxyB(i5,i2) = 1; 
        RxyB(i2,i1) = 1; 
        
    else %k is odd        
        
        if k==5 && m==r %last plaquette
            lmin = 2;
            i4 = ind(r,lmin-1);
        end
        
        H  (i1,i0) = jj(r0,r1);
        Rxx(i1,i0) = jj(r0,r1)*dd(r0,r1,1)^2;
        Rxy(i1,i0) = jj(r0,r1)*dd(r0,r1,1)*dd(r0,r1,2);
        
        H  (i2,i0) = jj(r0,r2);
        Rxx(i2,i0) = jj(r0,r2)*dd(r0,r2,1)^2;
        Rxy(i2,i0) = jj(r0,r2)*dd(r0,r2,1)*dd(r0,r2,2);
        
        H  (i2,i3) = jj(r3,r2);
        Rxx(i2,i3) = jj(r3,r2)*dd(r3,r2,1)^2;
        Rxy(i2,i3) = jj(r3,r2)*dd(r3,r2,1)*dd(r3,r2,2); 
        
        RxyB(i0,i4) = 1; 
        RxyB(i4,i3) = 1; 
        RxyB(i3,i0) = 1; 
        RxyA(i1,i5) = 1;
        RxyA(i5,i2) = 1; 
        RxyA(i2,i1) = 1; 
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
            Rxx(i4,i1) = jj(r4,r1)*dd(r4,r1,1)^2; 
            Rxy(i4,i1) = jj(r4,r1)*dd(r4,r1,1)*dd(r4,r1,2);
        else %k is odd
            H(  i1,i4) = jj(r4,r1);   
            Rxx(i1,i4) = jj(r4,r1)*dd(r4,r1,1)^2; 
            Rxy(i1,i4) = jj(r4,r1)*dd(r4,r1,1)*dd(r4,r1,2);
        end
    end 
end

H0 = H;
Rxx0 = Rxx;
Rxy0 = Rxy;
RxyA0 = RxyA;
RxyB0 = RxyB;

%GSen = 0;
GSflag = false;

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes
        
%% Put in fluxes
for pind = 1:numel(pv)
    
    for tind = 1:pv(pind)
    
        p = pind-1;
    
        H = H0;
        Rxx = Rxx0;
        Rxy = Rxy0;
        RxyA = RxyA0;
        RxyB = RxyB0;

        %If there's more fluxes asked for than plaquettes reduce to max
        if p > numplaq
            disp(['p too large, using max fluxes: ',num2str(numplaq)])
            p = numplaq;
        end
        %If we are setting more than half then change the base gauge
        if p > numplaq/2
            %set the base gauge to all fluxes
            [H,Rxx,Rxy,RxyA,RxyB] = flipAll(H,Rxx,Rxy,RxyA,RxyB);

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
            [H,Rxx,Rxy,RxyA,RxyB] = add_flux(ru,ku,mu, H,Rxx,Rxy,RxyA,RxyB);

        end

        %add_flux(1,0,1)

        RxyA = RxyA - RxyA';
        RxyB = RxyB - RxyB';

        %% Compute spectra

        %memory
       % tic
        [En,Ev,Ddd,~,Ixx,Ixy,Ixy2] = dos1_loop(Tv,H,Rxx,Rxy,RxyA,RxyB,bins,emax);
       % toc
        %comes out with a cell index corresponding to channels,
        % and in each, indices (Tbin,Ebin) 
        % We store it here with an additional flux config index 
        
        if nargin < 7 && ~GSflag %don't have the GS energy yet
            if p==0
                GSen = En;
                GSflag = true;
            else
                error('You need to have 0 flux first in the list or a given GSen')
            end
        end
        
 %       try
        bf = bfactor(p,En,Tv); %A vector of boltzmann factors along index Tbin
 %       catch
  %         bf = bfactor(p,En,Tv);
 %       end
        
        It{tind,1} = En;
        
        %Multiply along Tbin by the boltmann factors.
%try
% catch
 %    disp(size(Ddd)) 
% end
        It{tind,2} = bsxfun(@times,bf.', repmat(Ddd.',numel(Tv),1) );
        It{tind,3} = bsxfun(@times,bf.', Ixx  );   
        It{tind,4} = bsxfun(@times,bf.', Ixy  );
        It{tind,5} = bsxfun(@times,bf.', Ixy2 );  

        tmp(:,tind) = bf;
        Endist{pind} = [Endist{pind},En];
        
        if( any( isnan(bf) | (bf == Inf) | (bf == -Inf) ) )
            error('bf exploded!')
        end
        
    end
    
    % Reduce the data to something a bit more manageable
    if pv(pind) > 0
        toc
        disp([pind,pv(pind)])
        tic  
        
        for kk=1:knum
            Itt{pind,kk} = cell_mean( It(1:pv(pind),kk) );
%             if max(size(Itt{pind,kk})) == 0
%                 warning('confusion')
%             end
            Its{pind,kk} = cell_std ( It(1:pv(pind),kk) ) /sqrt(rp);
            if (kk>1) && any(any(Itt{pind,kk} < 0))
                error('wtf1')
            end
        end
    bfs(:,pind) = mean(tmp(:,1:pv(pind)),2); 
    bfserr(:,pind) = std(tmp(:,1:pv(pind)),1,2)./sqrt(pv(pind));
    %Ens(pind) = mean(Etmp(1:pv(pind)));
    %Enserr(pind) = std(Etmp(1:pv(pind)));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nested functions

    function mn = cell_mean(ce)
        sum = zeros(size(ce{1}));
        for in = 1:numel(ce)
            sum = sum + ce{in};
        end
        mn = sum/numel(ce);
    end

    function std = cell_std(ce)
        mn = cell_mean(ce);
        sum = zeros(size(ce{1}));
        for in = 1:numel(ce)
            sum = sum + (ce{in}-mn).^2;
        end
        std = sqrt(sum);
    end

    function bf = bfactor(p3,E3,T3)
        % lnbin = log( nn!/((nn-kk)!*kk!) ) = log( nchoosek(nn,kk) )
        lnbin = gammaln(numplaq+1) - gammaln(numplaq-p3+1) - gammaln(p3+1);
        bf = exp(lnbin-(E3-GSen)./T3 - pf1); 
    end

    function [H,Rxx,Rxy,RxyA,RxyB] = flipAll(H,Rxx,Rxy,RxyA,RxyB)
        
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
                    i5 = ind(r-1,l-1-2*k);
                    i7 = ind(r,l-2);
                    i8 = ind(r-1,l-3-2*k);
                    if k == 0 && m == 1
                        lmax = 6*(r^2 - (r-1)^2);
                        i3 = ind(r-1,lmax);
                        i8 = ind(r-1,lmax-1);
                        lmax = 6*((r+1)^2 - r^2);
                        i7 = ind(r,lmax);
                    end
                    if [r,k,m] == [1,1,1] %a random subcase
                        lmax = 6*(r^2 - (r-1)^2);
                        i8 = ind(r-1,lmax);
                    end

                    %Just need positions for drawing fluxes in the lattice pic
                    r0 = rr(r,k,m,0);     
     %               r1 = rr(r,k,m,1);     draw_text(r1,i1);
                    r2 = rr(r,k,m-1,1);   
                    r3 = rr(r-1,k,m-1,0); 
                   % r4 = rr(r-1,k,m,0);  draw_text(r0,i0);
                    r5 = rr(r-1,k,m-1,1);   
                    r7 = rr(r,k,m-1,0);   
                    r8 = rr(r-1,k,m-2,1);              

                    %Change the bond for every other plaquette 
                        %going around each ring level
                    sn  = -(-1)^m;
                    sn2 = -(-1)^r * sn;    %switch for odd/even k, but only if 
                                            % r is odd.
                    if mod(k,2) == 0
                        H  (i3,i2) = sn*H  (i3,i2);
                        Rxx(i3,i2) = sn*Rxx(i3,i2);
                        Rxy(i3,i2) = sn*Rxy(i3,i2);

                        RxyA(i3,i0) = sn*RxyA(i3,i0);
                        RxyB(i5,i2) = sn*RxyB(i5,i2);
                        RxyA(i7,i3) = sn*RxyA(i7,i3);
                        RxyB(i2,i8) = sn*RxyB(i2,i8);
                    else
                        H  (i2,i3) = -sn2*H  (i2,i3);
                        Rxx(i2,i3) = -sn2*Rxx(i2,i3);
                        Rxy(i2,i3) = -sn2*Rxy(i2,i3);

                        RxyB(i3,i0) = -sn2*RxyB(i3,i0);
                        RxyA(i5,i2) = -sn2*RxyA(i5,i2);
                        RxyB(i7,i3) = -sn2*RxyB(i7,i3);
                        RxyA(i2,i8) = -sn2*RxyA(i2,i8);                      
                    end

                end
            end
        end
        %set the center plaquette with the gauge string ending at the k=2 site
        %This amounts to changing the z-bond on the right side of the plaq.
        [H,Rxx,Rxy,RxyA,RxyB] = add_flux(0,0,1, H,Rxx,Rxy,RxyA,RxyB);
        
    end

    function [H,Rxx,Rxy,RxyA,RxyB] = add_flux(rrc,k,m, H,Rxx,Rxy,RxyA,RxyB)  
        %This function edits H and R
        if rrc==0
            [H,Rxx,Rxy,RxyA,RxyB] = add_flux(1,5,1, H,Rxx,Rxy,RxyA,RxyB);
            
            %Still need to change the top-left bond of 0-plaquette
            i0 = 3;
            i2 = 3;
            i1 = 1;
            i3 = 2;
            i6 = 4;
            i7 = 11;
            
            %These are only needed to draw the lattice picture!
            %Just need positions for drawing fluxes in the lattice pic
            r0 = rr(0,0,0,0);     
            r1 = rr(0,0,0,1);    
            r2 = rr(0,0,-1,1);   
            r3 = rr(-1,0,-1,0);   
            %r4 = rr(r-1,k,m,0);  
            %r5 = rr(r-1,k,m,1);   %odd
            rc = 0;
            if rmax>rc
                r6 = rr(1,0,0,1);     
                r7 = rr(0,0,-1,0);   
            end
            %r8 = rr(r-1,k,m-1,1);  
            
            %rc = 0;
            [H,Rxx,Rxy,RxyA,RxyB] = flip_0_2_bond(H,Rxx,Rxy,RxyA,RxyB);  
            
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
                r0 = rr(rc,k,m,0);     
                r1 = rr(rc,k,m,1);    
                r2 = rr(rc,k,m-1,1);   
                r3 = rr(rc-1,k,m-1,0);
                %r4 = rr(r-1,k,m,0);  draw_text(r4,i4)
                %r5 = rr(r-1,k,m,1);   draw_text(r5,i5)

                if rc < rmax
                    i6 = ind(rc+1,l+2*k+1);
                    i7 = ind(rc,l-2);
                    r6 = rr(rc+1,k,m,1);  
                    r7 = rr(rc,k,m-1,0);  
                    %r8 = rr(r-1,k,m-1,1);   
                end

                [H,Rxx,Rxy,RxyA,RxyB] = flip_0_2_bond(H,Rxx,Rxy,RxyA,RxyB);

            end    
        end
        
        function [H,Rxx,Rxy,RxyA,RxyB] = flip_0_2_bond(H,Rxx,Rxy,RxyA,RxyB)

            %Change everything originating from the (2,0) bond!
            %(0,2)
            if mod(k,2) == 0
                H  (i0,i2) = -H  (i0,i2);
                Rxx(i0,i2) = -Rxx(i0,i2);
                Rxy(i0,i2) = -Rxy(i0,i2); 
                RxyA(i3,i0) = -RxyA(i3,i0);
                RxyB(i2,i1) = -RxyB(i2,i1); 
                if rc < rmax
                    %disp([rc,rmax])
                    RxyA(i0,i7) = -RxyA(i0,i7);
                    RxyB(i6,i2) = -RxyB(i6,i2);   
                end
            else
                H  (i2,i0) = -H  (i2,i0);
                Rxx(i2,i0) = -Rxx(i2,i0);
                Rxy(i2,i0) = -Rxy(i2,i0);  

                RxyB(i3,i0) = -RxyB(i3,i0);
                RxyA(i2,i1) = -RxyA(i2,i1); 
                if rc < rmax
                    RxyB(i0,i7) = -RxyB(i0,i7);
                    RxyA(i6,i2) = -RxyA(i6,i2);   
                end
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

end