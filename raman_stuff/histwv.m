function [histw, histv] = histwv(v, w, min, max, bins)
    %Inputs: 
    %vv - values
    %ww - weights
    %minV - minimum value
    %maxV - max value
    %bins - number of bins (inclusive)
    
    %Outputs:
    %histw - wieghted histogram
    %histv (optional) - histogram of values    
   
    delta = (max-min)/(bins-1);
    subs = floor((v-min)/delta)+1;
    
    if nargout == 2
        histv = accumarray(subs(:),1,[bins,1]);
    end
    histw = accumarray(subs(:),w(:),[bins,1]);
end