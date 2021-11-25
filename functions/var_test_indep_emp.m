function [teststat, teststats_boot, A, c, shocks, H, H_cumul, C, C_cumul, allmins, loglik] = var_test_indep_emp(Y, Z, p, pml_settings, numboot, verbose, initSetting, initSettingBoot)

% Adapted from var_test_indep.m
% Replicates the approach of the empirical example in
% Gourieroux, Monfort, Renne (2017).

%% Initialization settings
if ~exist('initSetting', 'var')
    initSetting = 'cumul';  % Set default to fourth order cumulant
end

if ~exist('initSettingBoot', 'var')
    initSettingBoot = initSetting;  % Set default to the full-sample setting
end

%% Estimate SVAR

[T,n] = size(Y);

% Reduced form
[A, Sigma, c, res, cZ]  = var_estim2(Y, p, Z);

% Preliminary estimate of H from fourth order cumulants
[H_cumul,C_cumul] = ica_cumul(res);

if isa(initSetting, 'char')
    if strcmp(initSetting, 'cumul')
        [H, C, allmins, loglik] = pml(res, @(X) normal_mixture(X, pml_settings.mus, ...
            pml_settings.sigmas, pml_settings.ps), ...
            C_cumul, pml_settings.opts, 0);
        
    elseif strcmp(initSetting, 'GlobalSearch')
        [H, C, allmins, loglik] = pml(res, @(X) normal_mixture(X, pml_settings.mus, ...
            pml_settings.sigmas, pml_settings.ps), ...
            C_cumul, pml_settings.opts, 1);
        
    elseif strcmp(initSetting, 'MultiStart')
        [H, C, allmins, loglik] = pml(res, @(X) normal_mixture(X, pml_settings.mus, ...
            pml_settings.sigmas, pml_settings.ps), ...
            C_cumul, pml_settings.opts, 2);
        
    else
        error('Invalid char: Enter valid initialization setting...')
    end
    
elseif isa(initSetting, 'double')
    % Estimate H by PML. Input spherical residuals into pml routine.
    
    [H, C, allmins, loglik] = pml(res, @(X) normal_mixture(X, pml_settings.mus, ...
        pml_settings.sigmas, pml_settings.ps), ...
        initSetting, pml_settings.opts, 0);
else
    error('Invalid char: Enter valid initialization setting...')
end
% Shock estimates
shocks = res/(H');


%% Test statistic

corr_shocks_sq = corr(shocks.^2); % Correlation matrix for squared shocks
teststat       = sqrt((sum(corr_shocks_sq(:).^2)-n)/(n^2-n)); % Root mean squared off-diagonal correlation


if ~isempty(Z)  % Don't run bootstrap if there is an exogenous regressor in the VAR
    teststats_boot = NaN;  
else
    %%  Bootstrap
    
    teststats_boot = nan(numboot,1);
    
    if ~isempty(numboot) && numboot>0
        
        if verbose
            disp('Bootstrapping...');
        end
        
        rand_seeds = randi(2^32-1,numboot,1);
        timer = tic;
        
        if isa(initSettingBoot, 'char')
            if strcmp(initSettingBoot, 'MLE')
                initSettingBoot = C;
            end
        end
        
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
            teststats_boot(ib) = var_test_indep_emp(the_Y_boot, Z, p, pml_settings, 0, false, initSettingBoot);

            
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
end