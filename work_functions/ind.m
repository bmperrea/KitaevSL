function in = ind(r,l)
    in = 6*r^2 + l;
    %This puts 1 -> 1 , 2->1 , 3->2 , 4->2 ,...
    in = floor( (in+1)/2 );
end