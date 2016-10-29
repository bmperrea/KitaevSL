function [ps,rs,ks,ms] = get_NNNs(r,k,m,~)

    %The flag being 'p' says that you want ps insted of r,k,m s

    %The list of directions to make the NNN pattern
    dir = [1,3,4,5,0,1, 1,2,3,3,4,4,5,5,0,0,1,1 ,1];
    
    %initialize the ps
 %   if nargin > 3 && flag == 'p'
        ps = zeros(1,numel(dir));
 %   else
        rs = zeros(1,numel(dir));
        ks = zeros(1,numel(dir));
        ms = zeros(1,numel(dir));
 %   end

    for j=1:numel(dir)
        %store
%        if nargin > 3 && flag == 'p'
            ps(j) = 3*r*(r-1) + r*k + m + 1;
%        else
            rs(j) = r;
            ks(j) = k;
            ms(j) = m;
%        end
        %move to next
        [r,k,m] = move(r,k,m,dir(j));       
    end
    
%    if nargin > 3 && flag == 'p'
%        out = ps;
%    else
%        out = [rs;ks;ms];
%    end

end

function [r,k,m] = move(r,k,m,dir)
   %dir specifies the direction the same as k, 0 is up-right, 1 is right..
   
    %dr is the relative direction to the k-outward direction
    dr = dir - k;
    if dr<0
        dr = dr+6;
    end
    
    if r>0
    
        if dr == 0 %up-right
            r = r+1;
        elseif dr == 1 %right
         %   if m<r
                r=r+1;
                m=m+1;
          %  else %m==r
          %      r=r+1;
             %   k=k+1;
          %  end
        elseif dr == 2 %down-right
            if m<r
                m=m+1;
            else %m==r
                r=r+1;
                m=1;
                k=k+1; %TODO, check k==5
            end
        elseif dr == 3 %down-left
            if m<r
                r=r-1;
            else %m==r
                m=1;
                k=k+1; %TODO, check k==5
            end        
        elseif dr == 4 %left
            if m>1
                r = r-1;
                m = m-1;
            else %m==0
                r=r-1;
                m=r;
                k=k-1; %TODO, check k==0
            end        
        elseif dr == 5 %up-left
            if m>1
                m=m-1;
            else %m==1
                k=k-1; %TODO, check k==0
                m=r;
            end        
        end
              
    else %r==0
        
        r=1;
        m=1;
        k=dir-1;
        
    end
    
    %Watch overflow of k
    if k==6
        k=0;
    elseif k==-1;
        k=5;
    end  
    
end


%     Get l from r,k,m
%     l = (2*r+1)*k + 2*m;

%   Get r,k,m from l
%     ru = ceil( ( -1+sqrt(1+(plaq-1-eps)*4/3) )/2 );
%     iu = plaq - 1 - 3*ru*(ru-1);
%     ku = floor( (iu - 1 + eps)/ru );
%     mu = iu - ru*ku;