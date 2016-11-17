function [H,Rxx,Rxy,RxyA,RxyB,pbool] = set_gauge(rmax,p,H0,Rxx0,Rxy0,RxyA0,RxyB0)

    H = H0;
    Rxx = Rxx0;
    Rxy = Rxy0;
    RxyA = RxyA0;
    RxyB = RxyB0;

    %Generate a random flux pattern with p fluxes
    %ms = []; ks = ms; rs = ms;
    inds = [];
    numplaq = 3*rmax*(rmax+1)+1;
    while numel(inds) < p
        rans = randi(numplaq,p-numel(inds),1);
        inds = [inds;rans];
        %Remove duplicates
        inds = unique(inds);

    end
    pbool = false(1,numplaq);
    pbool(inds) = true;
    
    %Put in fluxes
    for nn = 1:p
    [ru,ku,mu] = get_polar(inds(nn));
    
    [H,Rxx,Rxy,RxyA,RxyB] = flip_plaquette(ru,ku,mu,rmax,...
                    H,Rxx,Rxy,RxyA,RxyB);              
    end

end

    
    
    