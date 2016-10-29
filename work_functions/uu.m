function u = uu(r1)
    vec = r1;
    x = vec(1); y = vec(2);
    u = [2*x*y,x^2-y^2];
end