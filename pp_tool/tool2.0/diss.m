function [ diss ] = diss( x, y )

    % global dissimilarity at time frame
    % input: microstate --> microstato da confrontare [1 x nch]
    %        potential --> potenziale misurato ad ogni time frame [1 x nch]
    X_hat = mean(x);
    Y_hat = mean(y);
    u = x - X_hat;
    v = y - Y_hat;
    gfp_microstate = std(u);
    gfp_potential = std(v);
    arg = ((u/gfp_microstate)-(v/gfp_potential)).^2;
%     arg = (u-v).^2;
    diss = sqrt((1/size(v,1)).*sum(arg));
end

