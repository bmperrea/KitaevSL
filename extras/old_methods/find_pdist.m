function [ps,Es,dist,dists] = find_pdist( rmax ,b,s , Tv ,number, ~)
%number = 100?
%rmax   = 15 to 50
% b = 10;
% s = .04
% Tv is a vector of temperatures

    numplaq = 1 + 3*rmax*(rmax+1);
    ps = round( linspace(0,numplaq,number) );
    
    Es = zeros(1,numel(ps));
    dist = Es;
    dists = zeros(numel(ps),numel(Tv));
    
    if ps(1) ~= 0
        error('must start with 0 flux')
    end
    for pind = 1:numel(ps)
        p = ps(pind);
        disp(pind/numel(ps))
        [~,En]=stretch_2D_6(rmax,b,s,p,0,-1);
        Es(pind) = En;
        if pind == 1
            GSen = En;
        end
        
        bf = bfactor(p,En,Tv);
        if( any( isnan(bf) | (bf == Inf) | (bf == -Inf) ) )
            error('bf exploded!')
        end
        dists(pind,:) = bf;
    end
    sums = sum( dists , 1 );
    dist = mean( bsxfun(@rdivide,dists,sums) , 2);
    
if nargin < 6 %If you provide a sixth flag argument we don't plot 
    hh=figure;
    plot(ps,Es);
    xlabel('p');
    ylabel('E/J');
    filename = ['stretch_Ensp_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);
    
    hh=figure;
    plot(ps,dist);
    xlabel('p');
    ylabel('dist');
    filename = ['stretch_distp_rmax_',num2str(rmax),'_b_',num2str(round(b)),...
        '_c_',num2str(round(1000*s))];
    saveas(hh,filename)
    print(hh, '-dpng', filename);
    print(hh, '-depsc', filename);
end
    
    function bf = bfactor(p3,E3,T3)
        % lnbin = log( nn!/((nn-kk)!*kk!) ) = log( nchoosek(nn,kk) )
        lnbin = gammaln(numplaq+1) - gammaln(numplaq-p3+1) - gammaln(p3+1);
        bf = exp(lnbin-(E3-GSen)./T3); 
    end

end


