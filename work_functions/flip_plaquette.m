function [H,Rxx,Rxy,RxyA,RxyB] = flip_plaquette(rrc,k,m,rmax,H,Rxx,Rxy,RxyA,RxyB)  
        %This function edits H and R
        if rrc==0
            [H,Rxx,Rxy,RxyA,RxyB] = flip_plaquette(1,5,1,rmax,H,Rxx,Rxy,RxyA,RxyB);
            
            %Still need to change the top-left bond of 0-plaquette
            i0 = 3;
            i2 = 3;
            i1 = 1;
            i3 = 2;
            i6 = 4;
            i7 = 11;
            
            %These are only needed to draw the lattice picture!
            %Just need positions for drawing fluxes in the lattice pic
%             r0 = rr(0,0,0,0);     %draw_text(r0,i0)
%             r1 = rr(0,0,0,1);    % draw_text(r1,i1)
%             r2 = rr(0,0,-1,1);   % draw_text(r2,i2)
%             r3 = rr(-1,0,-1,0);  % draw_text(r3,i3)
            %r4 = rr(r-1,k,m,0);  
            %r5 = rr(r-1,k,m,1);   %odd
            rc = 0;
%             if rmax>rc
%                 r6 = rr(1,0,0,1);   %  draw_text(r6,i6)
%                 r7 = rr(0,0,-1,0);  %  draw_text(r7,i7)
%             end
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
            %    r0 = rr(rc,k,m,0);    % draw_text(r0,i0)
             %   r1 = rr(rc,k,m,1);    % draw_text(r1,i1)
            %    r2 = rr(rc,k,m-1,1);  % draw_text(r2,i2)
             %   r3 = rr(rc-1,k,m-1,0); %draw_text(r3,i3)
                %r4 = rr(r-1,k,m,0);  draw_text(r4,i4)
                %r5 = rr(r-1,k,m,1);   draw_text(r5,i5)

                if rc < rmax
                    i6 = ind(rc+1,l+2*k+1);
                    i7 = ind(rc,l-2);
            %        r6 = rr(rc+1,k,m,1);  % draw_text(r6,i6)
            %        r7 = rr(rc,k,m-1,0); %  draw_text(r7,i7)
                    %r8 = rr(r-1,k,m-1,1);   
                end

                flip_0_2_bond()

            end    
        end
        
        function flip_0_2_bond()

            %Change everything originating from the (2,0) bond!
            %(0,2)
            if mod(k,2) == 0
                H  (i0,i2) = -H  (i0,i2);
                Rxx(i0,i2) = -Rxx(i0,i2);
                Rxy(i0,i2) = -Rxy(i0,i2);  
               % draw_line(r0,r2,H(i0,i2));
 %          try
                RxyA(i3,i0) = -RxyA(i3,i0);
%           catch
 %              disp('yaya')
%           end
               % draw_line2(r3,r0,RxyA(i3,i0))
                RxyB(i2,i1) = -RxyB(i2,i1); 
              %  draw_line2(r2,r1,RxyB(i2,i1))
                if rc < rmax
                    %disp([rc,rmax])
                    RxyA(i0,i7) = -RxyA(i0,i7);
                  %  draw_line2(r0,r7,RxyA(i0,i7))
                    RxyB(i6,i2) = -RxyB(i6,i2); 
                   % draw_line2(r6,r2,RxyB(i6,i2))     
                end
            else
                H  (i2,i0) = -H  (i2,i0);
                Rxx(i2,i0) = -Rxx(i2,i0);
                Rxy(i2,i0) = -Rxy(i2,i0);  
              %  draw_line(r2,r0,H(i2,i0)); 

                RxyB(i3,i0) = -RxyB(i3,i0);
               % draw_line2(r3,r0,RxyB(i3,i0))
                RxyA(i2,i1) = -RxyA(i2,i1); 
               % draw_line2(r2,r1,RxyA(i2,i1))
                if rc < rmax
                    RxyB(i0,i7) = -RxyB(i0,i7);
                   % draw_line2(r0,r7,RxyB(i0,i7))
                    RxyA(i6,i2) = -RxyA(i6,i2); 
                   % draw_line2(r6,r2,RxyA(i6,i2))   
                end
            end
        
        end
        
    end