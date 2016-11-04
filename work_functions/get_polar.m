function [ru,ku,mu,ju,iu] = get_polar(plaq)
    %Convert index number to polar coordinate data
    eps = 0.0001; %A small number to help avoid numerical errors with ceil
    ru = ceil( ( -1+sqrt(1+(plaq-1-eps)*4/3) )/2 );
    iu = plaq - 1 - 3*ru*(ru-1);
    ku = floor( (iu - 1 + eps)/ru );
    mu = iu - ru*ku;
    if ru==0
        mu = 0;
        ku = 0;
    end
    ju = ru*(ru-1)/2 + mu + 1;
end