function [A, Sigma, c, res, cZ, aic, sc, hq] = var_estim(Y, p, Z)
% Modified 6/15/21. Add argument for exogenous variables Z.

    %% VAR least-squares estimation
    
    nZ = size(Z, 2);      % Number of exogenous variables 
    k  = size(Y, 2);      % Number of endogenous variables
    T  = size(Y, 1) - p;  % Estimated sample size 
    
    X     = [lagmatrix(Y, 1:p) Z ones(size(Y,1),1)];  
    beta  = (X(p+1:end,:)\Y(p+1:end,:))';
    A     = beta(:,1:end-1-nZ);
    c     = beta(:,end);
    cZ    = beta(:, end-nZ:end-1);            % Coefficient on endogenous regressor
    res   = Y(p+1:end,:)-X(p+1:end,:)*beta';
    Sigma = (res'*res)/(T- k*p - nZ - 1);  % Degrees of freedom adjustment
    
    
    %% Compute Akaike Information Criterion (aic), Schwarz Criterion (sc),
    % and Hannan-Quinn Criterion (hq)

    l       = -T/2 * (k * (1+log(2*pi)) + log(det(cov(res, 1)))); % log-likelihood
    nParams = k*(p*k + 1 + nZ);                     % Number of estimated parameters
    
    % See EViews documentation for details
    aic = -2*l/T + 2 *nParams/ T;
    sc  = -2*l/T + nParams*log(T) / T;
    hq  = -2*l/T + 2 * nParams * log(log(T)) /T;

end