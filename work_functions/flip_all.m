function [H,Rxx,Rxy,RxyA,RxyB]=flip_all(rmax,H,Rxx,Rxy,RxyA,RxyB)
       %set the base gauge to all fluxes
    for r=1:rmax
        for k=0:5
            for m = 1:r
                %Change the (2,3) bond and the associated NNN bonds
                
                l = (2*r+1)*k + 2*m;

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
               % r0 = rr(r,k,m,0);     %draw_text(r0,i0);
 %               r1 = rr(r,k,m,1);     draw_text(r1,i1);
              %  r2 = rr(r,k,m-1,1);   %draw_text(r2,i2);
              %  r3 = rr(r-1,k,m-1,0); %draw_text(r3,i3);
               % r4 = rr(r-1,k,m,0);  draw_text(r0,i0);
              %  r5 = rr(r-1,k,m-1,1);  % draw_text(r5,i5);
              %  r7 = rr(r,k,m-1,0);  % draw_text(r7,i7);
              %  r8 = rr(r-1,k,m-2,1); %draw_text(r8,i8);              

                %Change the bond for every other plaquette 
                    %going around each ring level
                sn  = -(-1)^m;
                sn2 = -(-1)^r * sn;    %switch for odd/even k, but only if 
                                        % r is odd.
                if mod(k,2) == 0
                    H  (i3,i2) = sn*H  (i3,i2);
                    Rxx(i3,i2) = sn*Rxx(i3,i2);
                    Rxy(i3,i2) = sn*Rxy(i3,i2);
                 %   draw_line(r3,r2,H(i3,i2))
                  
                    RxyA(i3,i0) = sn*RxyA(i3,i0);
                   % draw_line2(r3,r0,RxyA(i3,i0))
                    RxyB(i5,i2) = sn*RxyB(i5,i2);
                   % draw_line2(r5,r2,RxyB(i5,i2))
                    RxyA(i7,i3) = sn*RxyA(i7,i3);
                  %  draw_line2(r7,r3,RxyA(i7,i3))
                    RxyB(i2,i8) = sn*RxyB(i2,i8);
                  %  draw_line2(r2,r8,RxyB(i2,i8))
                else
                    H  (i2,i3) = -sn2*H  (i2,i3);
                    Rxx(i2,i3) = -sn2*Rxx(i2,i3);
                    Rxy(i2,i3) = -sn2*Rxy(i2,i3);
                  %  draw_line(r2,r3,H(i2,i3))

                    RxyB(i3,i0) = -sn2*RxyB(i3,i0);
                  %  draw_line2(r3,r0,RxyB(i3,i0))
                    RxyA(i5,i2) = -sn2*RxyA(i5,i2);
                  %  draw_line2(r5,r2,RxyA(i5,i2))
                    RxyB(i7,i3) = -sn2*RxyB(i7,i3);
                  %  draw_line2(r7,r3,RxyB(i7,i3))
                    RxyA(i2,i8) = -sn2*RxyA(i2,i8); 
                  %  draw_line2(r2,r8,RxyA(i2,i8))                       
                end
                              
            end
        end
    end
    %set the center plaquette with the gauge string ending at the k=2 site
    %This amounts to changing the z-bond on the right side of the plaq.
    [H,Rxx,Rxy,RxyA,RxyB] = flip_plaquette(0,0,0,rmax,H,Rxx,Rxy,RxyA,RxyB);
end