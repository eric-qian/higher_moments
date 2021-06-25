function [teststat, teststats_boot, A, c, shocks, H, H_cumul, C, C_cumul] = var_test_indep_emp(Y, Z, p, pml_settings, numboot, verbose, findGlobal)

    % Adapted from var_test_indep.m
    % Replicates the approach of the empirical example in 
    % Gourieroux, Monfort, Renne (2017).
    
    %% Estimate SVAR
    
    [T,n] = size(Y);
        
    % Reduced form
    [A, Sigma, c, res, cZ]  = var_estim2(Y, p, Z);
    
    % Preliminary estimate of H from fourth order cumulants
    [H_cumul,C_cumul] = ica_cumul(res);

    % Estimate H by PML. Input spherical residuals into pml routine.
    [H, C] = pml(res, @(X) normal_mixture(X, pml_settings.mus, pml_settings.sigmas, pml_settings.ps), C_cumul, pml_settings.opts, findGlobal);
    
    % Shock estimates
    shocks = res/(H');

    
    %% Test statistic
    
    corr_shocks_sq = corr(shocks.^2); % Correlation matrix for squared shocks
    teststat       = sqrt((sum(corr_shocks_sq(:).^2)-n)/(n^2-n)); % Root mean squared off-diagonal correlation
    teststats_boot = NaN;  % Bootstrap procedure is incomplete

end