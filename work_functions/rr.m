function r1 = rr(r,k,m,l1)
    interc = { [0,0]          , [sqrt(3),-1]/2, [sqrt(3),-3]/2 , [0,-2]         , [-sqrt(3),-3]/2 , [-sqrt(3),-1]/2 };
    sloper = { [-sqrt(3),3]/2 , [sqrt(3),3]/2 , [sqrt(3),0]    , [sqrt(3),-3]/2 , [-sqrt(3),-3]/2 , [-sqrt(3),0]    };
    slopem = { [sqrt(3),0]    , [sqrt(3),-3]/2, [-sqrt(3),-3]/2, [-sqrt(3),0]   , [-sqrt(3),3]/2  , [sqrt(3),3]/2   };
    slopel = { [sqrt(3),-1]/2 , [0,-1]        , [-sqrt(3),-1]/2, [-sqrt(3),1]/2 , [0,1]           , [sqrt(3),1]/2   };
    r1 = interc{k+1} + r*sloper{k+1} + m*slopem{k+1} + (l1)*slopel{k+1} + [0,1];
end