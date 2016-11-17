function [ru,ku,mu,ju,iu] = get_polar(plaq)
    %Convert index number to polar coordinate data
    %Input can be vector or matrix of plaq's and output will be
    %element-wise vectors or matrices of polor coord data
    eps = 0.0001; %A small number to help avoid numerical errors with ceil
    ru = ceil( ( -1+sqrt(1+(plaq-1-eps)*4/3) )/2 );
    iu = plaq - 1 - 3*ru.*(ru-1);
    ku = floor( (iu - 1 + eps)./ru );
    mu = iu - ru.*ku;
    zer = ru==0;
    mu(zer) = 0;
    ku(zer) = 0;
%     if ru==0
%         mu = 0;
%         ku = 0;
%     end
    ju = ru.*(ru-1)/2 + mu + 1;
end