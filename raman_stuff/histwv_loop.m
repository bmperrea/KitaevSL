function [histw, histv, subs, inds] = histwv_loop(v, w, bmin, bmax, bins, subs, inds)

    %Inputs: 
    %vv - values as a vector
    %ww - weights as a matrix
    %minV - minimum value
    %maxV - max value
    %bins - number of bins (inclusive)
    
    %Outputs:
    %histw - wieghted histogram
    %histv (optional) - histogram of values    
    if min(size(w)) == 1
        
        ins = (bmin<=v & v<bmax);
        if nargout == 1
        	histw = histwv(v(ins), w(ins), bmin, bmax, bins);
        else
            [histw,histv] = histwv(v(ins), w(ins), bmin, bmax, bins);
        end
        subs = -1;
        inds = -1;
        
%     elseif min(size(w)) == 2
%         
%         ins = (bmin<=v & v<bmax);
%         for in = 1:(numel(w)/numel(v))
%             if nargout == 1
%                 if size(v,1) == 1
%                     histw = histwv(v(ins), w(in,ins), bmin, bmax, bins);
%                 else
%                     histw = histwv(v(ins), w(ins,in), bmin, bmax, bins);
%                 end
%             else
%                 if size(v,1) == 1
%                     [histw,histv] = histwv(v(ins), w(in,ins), bmin, bmax, bins);
%                 else
%                     [histw,histv] = histwv(v(ins), w(ins,in), bmin, bmax, bins);
%                 end
%             end
%         end
%         subs = -1;
%         inds = -1;
        
    elseif size(v,1) == 1  %This is for multiple dimension for w...
    
        if nargin <= 5 || numel(subs) <= 1  || numel(inds) <= 1
            delta = (bmax-bmin)/(bins-1); %bin width
            subs = floor((v-bmin)/delta)+1; %bin number
            [subs,inds] = sort(subs); %Sorted bin numbers

            %Remove elements out of bounds and sort w and v (along v)
            spots = 1<=subs & subs<bins; %Could do this after sort...
            subs = subs( spots ); %could use inds directly...
            inds = inds( spots );
        end            
            
        w = w(:,inds);

        %Take the cumulative sum of w
        w = cumsum(w,2);

        %Find where subs changes to decide where to find the appropriate sum
        changes = [find(diff(subs)), numel(subs)]; %spot with last same value

        %Find the binning index at these changes
        subIns = subs(changes);

        %find the weights at these changes
        wC = diff( [zeros(size(w,1),1), w(:,changes)] ,1,2);

        %Set histw at these indices
        histw = zeros(size(w,1),bins);
        histw(:,subIns) = wC;
        
    elseif size(v,2) == 1 
        
        if nargin <= 5 || numel(subs) <= 1  || numel(inds) <= 1
            delta = (bmax-bmin)/(bins-1); %bin width
            subs = floor((v-bmin)/delta)+1; %bin number
            [subs,inds] = sort(subs); %Sorted bin numbers

            %Remove elements out of bounds and sort w and v (along v)
            spots = 1<=subs & subs<bins; %Could do this after sort...
            subs = subs( spots ); %could use inds directly...
            inds = inds( spots );
        end
        
        w = w(inds,:);

        %Take the cumulative sum of w
        w = cumsum(w,1);

        %Find where subs changes to decide where to find the appropriate sum
        changes = [find(diff(subs)); numel(subs)]; %spot with last same value

        %Find the binning index at these changes
        subIns = subs(changes);

        %find the weights at these changes
        wC = diff( [zeros(1,size(w,2)); w(changes,:)] ,1,1);

        %Set histw at these indices
        histw = zeros(size(w,2),bins);
        histw(:,subIns) = wC.';
        
    else
        error('Something is probably wrong with the input');
    end
        
    if nargout > 1 && ~exist('histv','var')
        if ~exist('subs','var')
            delta = (bmax-bmin)/(bins-1); %bin width
            subs = floor((v-bmin)/delta)+1; %bin number
        end
        histv = accumarray(subs(:),1,[bins,1]);
    end
    
end