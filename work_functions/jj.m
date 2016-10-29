function j = jj(r1,r2,b)
    j = 1 + b - b*norm( dd(r1,r2) );      
    if j <= 0
        error('j is gone!')
    end
end