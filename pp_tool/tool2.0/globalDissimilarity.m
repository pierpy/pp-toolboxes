function diss = globalDissimilarity(X, Y)
    % global dissimilarity at time frame
    % input: microstate --> microstato da confrontare [1 x nch]
    %        potential --> potenziale misurato ad ogni time frame [1 x nch]
    X_hat = mean(X);
    Y_hat = mean(Y);
    u = X - X_hat;
    v = Y - Y_hat;
    gfp_microstate = std(u);
    gfp_potential = std(v);
    arg = ((u/gfp_microstate)-(v/gfp_potential)).^2;
diss = sqrt((1/size(v,1)).*sum(arg));

