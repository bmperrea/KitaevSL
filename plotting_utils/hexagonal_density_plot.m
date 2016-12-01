function hexagonal_density_plot(rmax,s,phist)
    
    c = s/(sqrt(3)*rmax+sqrt(3)/2);
    numplaq = 3 * rmax * (rmax +1 ) +1;
    
%   Here's the labeling of sites:
%             -     -
%          -     6     -
%          7     0     -
%             2     1     -
%             3     4     -
%          8     5     -

    %figure;
    axis tight
    colormap hot
    %axis square

    [rus,kus,mus] = get_polar(1:numplaq);
    
    pdensity = reshape( phist,numplaq,1);% / max(phist) 
    
    Xs = zeros(numplaq,6);
    Ys = zeros(numplaq,6);
    
    for jj= 1:numplaq

        r = rus(jj);
        k = kus(jj);
        m = mus(jj);    
        %l = (2*r+1)*k + 2*m;

        %The actual positions after strain
        r0 = rr(r,k,m,0);     %even sublattice
        r1 = rr(r,k,m,1);     %odd
        r2 = rr(r,k,m-1,1);   %odd
        r3 = rr(r-1,k,m-1,0); %even
        r4 = rr(r-1,k,m,0);   %even
        r5 = rr(r-1,k,m-1,1);   %odd

        vecs = [r0;r1;r4;r5;r3;r2];
        
        Xs(jj,:) = vecs(:,1).';
        Ys(jj,:) = vecs(:,2).';

    end
    
    patch(Xs.', Ys.', pdensity.')   
    colorbar
    
    function r1 = rr(r,k,m,l1)
        interc = { [0,0]          , [sqrt(3),-1]/2, [sqrt(3),-3]/2 , [0,-2]         , [-sqrt(3),-3]/2 , [-sqrt(3),-1]/2 };
        sloper = { [-sqrt(3),3]/2 , [sqrt(3),3]/2 , [sqrt(3),0]    , [sqrt(3),-3]/2 , [-sqrt(3),-3]/2 , [-sqrt(3),0]    };
        slopem = { [sqrt(3),0]    , [sqrt(3),-3]/2, [-sqrt(3),-3]/2, [-sqrt(3),0]   , [-sqrt(3),3]/2  , [sqrt(3),3]/2   };
        slopel = { [sqrt(3),-1]/2 , [0,-1]        , [-sqrt(3),-1]/2, [-sqrt(3),1]/2 , [0,1]           , [sqrt(3),1]/2   };
        r1 = interc{k+1} + r*sloper{k+1} + m*slopem{k+1} + (l1)*slopel{k+1} + [0,1];
        r1 = r1 + c*uu(r1);
    end

    function u = uu(r1)
        vec = r1;
        x = vec(1); y = vec(2);
        u = [2*x*y,x^2-y^2];
    end

end