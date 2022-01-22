function [teststat, teststats_boot, A, c, shocks, H, H_cumul, C, C_cumul, allmins, loglik] = var_test_indep(Y, Z, p, varargin)

% Based on the empirical example Gourieroux, Monfort, Renne (2017).
%
% Input:
% - Y:               Endogenous variables.
% - Z:               Exogenous variables.
% - p:               Number of lags.
% - pml_settings:    PML mixture settings.
% - numboot:         Number of bootstrap draws.
% - verbose:         =1 to show output progress, =0 otherwise.
% - initSetting:     Initalize using matrix, fourth-order cumulant
%                    "cumul") or Global Search algorithm
%                    "GlobalSearch").
% - initSettingBoot: Initialization setting for bootstrap draws. Same  
%                    options as initSetting. In addition, can input "full"
%                    to initialize at the full-sample value.
%
% Output:
% - teststat:        Test statistic.      
% - teststats_boot:  Test statistic for bootstrap draws.
% - A:               VAR coefficients.
% - c:               VAR constant.
% - shocks:          Shock estimates.
% - H:               Unwhitened mixing matrix.
% - H_cumul:         Unwhitened mixing matrix (fourth-order cumulant).
% - C:               Whitened mixing matrix.
% - C_cumul:         Whitened mixing matrix (fourth-order cumulant).
% - allmins:         If global optimization routine is run, stores ALL
%                    values of C.
% - loglik:          Log-liklihood.

%% Parse input

% Default settings
numboot_default         = 0;
verbose_default         = 0;
initSetting_default     = 'GlobalSearch';
initSettingBoot_default = 'full';

% PML default settings
pml_settings_default.ps         = repmat(0.5,2,3);                       % Mixture probabilities
pml_settings_default.mus        = [repmat( 0.1,1,3); repmat(-0.1,1,3)];  % Mixture means
pml_settings_default.sigmas     = [0.5 1.3 0.7; 1.32 0.54 1.22];         % Mixture std dev's
pml_settings_default.opts       = optimoptions('fmincon', 'Display', 'notify', ...
    'SpecifyObjectiveGradient',  true, ...
    'SpecifyConstraintGradient', true, ...
    'CheckGradients',            false, ...
    'MaxFunctionEvaluations',    5e3);                       % fmincon options


% Setup parser
parser = inputParser;
addRequired(parser, 'Y', @isnumeric);
addRequired(parser, 'Z', @isnumeric);
addRequired(parser, 'p', @isscalar);
addOptional(parser, 'pml_settings',    pml_settings_default,    @isstruct);
addOptional(parser, 'numboot',         numboot_default,         @isnumeric);
addParameter(parser,'verbose',         verbose_default,         @isscalar);
addParameter(parser,'initSetting',     initSetting_default);
addParameter(parser,'initSettingBoot', initSettingBoot_default);

% Unpack parameters
parse(parser,Y, Z, p, varargin{:});
pml_settings    = parser.Results.pml_settings;
numboot         = parser.Results.numboot;
verbose         = parser.Results.verbose;
initSetting     = parser.Results.initSetting;
initSettingBoot = parser.Results.initSettingBoot;

%% Estimate SVAR

[T,n] = size(Y);

% Reduced form
[A, ~, c, res, ~]  = var_estim(Y, p, Z);

% Preliminary estimate of H from fourth order cumulants
[H_cumul,C_cumul] = ica_cumul(res);

% Log density
f = @(X) normal_mixture(X, pml_settings.mus, ...
    pml_settings.sigmas, pml_settings.ps);


% Run PML routine
if isa(initSetting, 'char')  % Fourth order cumulant or Global Search
    if strcmp(initSetting, 'cumul')
        [H, C, allmins, loglik] = pml(res, f, C_cumul, pml_settings.opts, 'off');
        
        
    elseif strcmp(initSetting, 'GlobalSearch')
        [H, C, allmins, loglik] = pml(res, f, C_cumul, pml_settings.opts, initSetting);
        
    else
        error('Invalid char: Enter valid initialization setting...')
    end
    
elseif isa(initSetting, 'double')  % Initialize at given matrix
    % Estimate H by PML. Input spherical residuals into pml routine.
    
    [H, C, allmins, loglik] = pml(res, f, initSetting, pml_settings.opts, 'off');
else
    error('Invalid char: Enter valid initialization setting...')
end


% Shock estimates
shocks = res/(H');


%% Compute test statistic

corr_shocks_sq = corr(shocks.^2);                             % Correlation matrix for squared shocks
teststat       = sqrt((sum(corr_shocks_sq(:).^2)-n)/(n^2-n)); % Root mean squared off-diagonal correlation


%%  Bootstrap

if ~isempty(Z)  % Don't run bootstrap if there is an exogenous regressor in the VAR
    teststats_boot = NaN;
else
    
    teststats_boot = nan(numboot,1);
    
    if ~isempty(numboot) && numboot>0  % Run bootstrap
        
        if verbose
            disp('Bootstrapping...');
        end
        
        rand_seeds = randi(2^32-1,numboot,1);
        timer = tic;
        
        if isa(initSettingBoot, 'char')  % MLE initialization setting
            if strcmp(initSettingBoot, 'full')
                initSettingBoot = C;
            end
        end
        
        for ib=1:numboot
                        
            rng(rand_seeds(ib), 'twister'); % Seed RNG
            
            % Random block of initial values
            the_ind_init = randi(T-p);
            the_Y_init   = Y(the_ind_init:the_ind_init+p-1,:);
            
            
            % Draw bootstrap shocks *independently across columns*
            the_shocks_boot = shocks;
            for j=1:n
                the_shocks_boot(:,j) = the_shocks_boot(randi(T-p,T-p,1),j);
            end
            
            % Generate bootstrap sample
            the_Y_boot = var_sim(A, c, the_shocks_boot*H', T, the_Y_init);
            
            % Calculate test statistic on bootstrap sample
            teststats_boot(ib) = var_test_indep(the_Y_boot, Z, p, pml_settings, 0,...
                'verbose', false, 'initSetting', initSettingBoot);
            
            
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