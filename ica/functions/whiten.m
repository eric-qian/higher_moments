function [U, Sigma_chol] = whiten(X)

    % Data whitening
    
    Sigma = cov(X);
    Sigma_chol = chol(Sigma, 'lower');
    U = X/(Sigma_chol');

end