function [teststat, teststats_boot, A, c, shocks, H, H_cumul] = var_test_indep(Y, p, pml_settings, numboot, verbose)

    % Test independence of structural shocks estimated by PML
    % using residual-based VAR bootstrap that imposes the null hypothesis
    
    % MPM 2021-05-27

    
    %% Estimate SVAR
    
    [T,n] = size(Y);
    
    % Reduced form
    [A, ~, c, res] = var_estim(Y, p);
    
    % Preliminary estimate of H from fourth order cumulants
    [H_cumul,C_cumul] = ica_cumul(res);

    % Estimate H by PML
    H = pml(res, @(X) normal_mixture(X, pml_settings.mus, pml_settings.sigmas, pml_settings.ps), C_cumul, pml_settings.opts);
    
    % Shock estimates
    shocks = res/(H');

    
    %% Test statistic
    
    corr_shocks_sq = corr(shocks.^2); % Correlation matrix for squared shocks
    teststat = sqrt((sum(corr_shocks_sq(:).^2)-n)/(n^2-n)); % Root mean squared off-diagonal correlation
    
    
    %%  Bootstrap
    
    teststats_boot = nan(numboot,1);
    
    if ~isempty(numboot) && numboot>0
        
        if verbose
            disp('Bootstrapping...');
        end
        
        rand_seeds = randi(2^32-1,numboot,1);        
        timer = tic;
        
        for ib=1:numboot
            
%             clc
%             disp(['Running ' num2str(ib) ' of ' num2str(numboot) '...' ])
            
            rng(rand_seeds(ib), 'twister'); % Seed RNG

            % Random block of initial values
            the_ind_init = randi(T-p);
            the_Y_init = Y(the_ind_init:the_ind_init+p-1,:);

            % Draw bootstrap shocks *independently across columns*
            the_shocks_boot = shocks;
            for j=1:n
                the_shocks_boot(:,j) = the_shocks_boot(randi(T-p,T-p,1),j);
            end

            % Generate bootstrap sample
            the_Y_boot = var_sim(A, c, the_shocks_boot*H', T, the_Y_init);

            % Calculate test statistic on bootstrap sample
            teststats_boot(ib) = var_test_indep(the_Y_boot, p, pml_settings, 0, false);

            % Print progress
            if verbose && mod(ib,ceil(numboot/100))==0
                fprintf('%s%3d%s\n', repmat(' ',1,round(50*ib/numboot)), round(100*ib/numboot), '%');
            end

        end
        
        elapsed_time = toc(timer);
        
        if verbose
            disp('Done. Elapsed time (min):');
            disp(elapsed_time/60);
        end
        
    end

end