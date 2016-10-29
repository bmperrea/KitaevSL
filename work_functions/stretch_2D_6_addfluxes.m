function [En,H,Es] = stretch_2D_6_addfluxes(H,rmax,b,s,p,~,~)
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
%S = imax/2;


%% Put in fluxes

numplaq = 1 + 3*rmax*(rmax+1); %The number of plaquettes
%If there's more fluxes asked for than plaquettes reduce to max
if p > numplaq
    disp(['p too large, using max fluxes: ',num2str(numplaq)])
    p = numplaq;
end
%If we are setting more than half then change the base gauge
flip_all()

%Generate a random flux pattern with p fluxes
generate_pattern()
% if numel(ls)<p
%     error(['ls_',num2str(numel(ls))),'_p_',num2tr(p)])
% end

%Put in fluxes
put_in_fluxes()

Es = svd(H);
En = -sum(Es);


%% Nested functions

    function put_in_fluxes()
        
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
        %         l = 2*m + k*(2*r1+
         %   end
        end
        
    end

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

    function flip_all()
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
    end

    function generate_pattern()
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