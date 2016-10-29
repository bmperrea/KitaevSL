function inds = generate_pattern(rmax,p)
    %ms = []; ks = ms; rs = ms;
    inds = [];
    plaqs = 3*rmax*(rmax+1)+1;
    while numel(inds) < p
        rans = randi(plaqs,p-numel(inds),1);
        inds = [inds;rans];
        %Remove duplicates
        inds = unique(inds);
    end
end