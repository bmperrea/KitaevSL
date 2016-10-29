function d = dd(r1,r2,in,c)
%The third argument is optional
    du = ( uu(r1) - uu(r2) );
    d1 = r1 - r2 + c*du;
    if nargin == 3
        d = d1(in);
    else
        d = d1;
    end        
    if c*norm(du) > 1
        error('too much strain')
    end
end