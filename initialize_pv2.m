%Need to sample from a distribution representing the energy penalties
%  ~ exp(-E_p) as a function of p
numplaq = 1 + 3*rmax*(rmax+1);
[ps,Es,dist,dists] = find_pdist( rmax ,b,s , Tv , ceil(pnum/2.5));

%a rough estimate of the log of the partition function for each temp.
pf1 = log( sum(dists,1) );

tic

dist2 = interp1(ps,dist,0:numplaq);

subs = gendist(dist2,1,pnum);
pv = accumarray(subs.',1,[numplaq+1,1]);

if pv(1) == 0
    pv(1) = 1;
end
if pv(numplaq+1) == 0
    pv(numplaq+1) = 1;
end

%Bump up the value to at least 3 everywhere that we compute to make the
%error analysis easier/more reliable
pv(0<pv & pv<rp) = rp;

toc

figure;
plot(pv)