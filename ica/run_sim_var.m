clear all;
addpath('../functions');

% Simulate bootstrap test of mutual shock independence

% DGP:
% Y_t = c + A_1*Y_{t-1} + ... + A_p*Y_{t-p} + H*shocks_t,
% where:
%   shocks_{i,t} = sigma_t * zeta_{i,t},
%   sigma_t ~ log_normal (independent across t),
%   zeta_{i,t} ~ t(df_i) (independent across i and t),
%   sigma and zeta are independent of each other

% MPM 2021-05-26


%% Settings

% DGP
dfs    = [9 11 13];                         % df of zeta t-distributions
A_true = [0.8 0.3 0; 0 0.7 0; 0.2 0.2 0.2]; % True VAR(p) matrix [A_1,A_2,...,A_p]
H_true = [1 0 0; ones(1,3); 3 2 1];         % True impact impulse responses (diagonal elements should be positive)
tau    = 0;                                 % Std dev of log(sigma_t) (mean is chosen so that E[sigma_t^2]=1)
T      = 250;                               % Sample size
burnin = 100;                               % Simulation burn-in

% Estimation/test
p_estim = 1;    % Estimation lag length
numboot = 5e2;  % No. of bootstrap iterations

% PML mixture parameters
pml_settings.ps     = repmat(0.5,2,3);                       % Mixture probabilities
pml_settings.mus    = [repmat( 0.1,1,3); 
                       repmat(-0.1,1,3)];  % Mixture means
pml_settings.sigmas = [0.5 1.3 0.7; 
                       1.32 0.54 1.22];     % Mixture std dev's
pml_settings.opts = optimoptions('fmincon', 'Display', 'notify', ...
                        'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                        'CheckGradients', false, 'MaxFunctionEvaluations', 5e3); % fmincon options

% Simulation
quants = [0.68 0.75 0.9 0.95];  % Bootstrap quantiles to store
numsim = 2e2;                   % Number of Monte Carlo repetitions
rng(202105272, 'twister');      % Seed RNG


%% DGP

[n,np] = size(A_true);
p = np/n;

comp_form = [A_true; eye(np-n) zeros(np-n,n)]; % Companion form
comp_eigs = eig(comp_form);

disp('Eigenvalues of companion matrix (magnitudes)');
disp(sort(abs(comp_eigs),'descend'));


%% Simulate

poolobj = parpool('local', 20);

teststats             = nan(numsim,1);
teststats_boot_quants = nan(numsim,length(quants));
H_estims = nan(numsim,n,n,2); % Last dimension: 1 is PML and 2 is cumulant estimate
rand_seeds = randi(2^32-1,numsim,1);
timer = tic;

disp('Simulating...');

parfor ii=1:numsim
    
    rng(rand_seeds(ii), 'twister'); % Seed RNG
    
    % Simulate VAR data
    zeta   = trnd(repmat(dfs,T+burnin-p,1))./sqrt(dfs./(dfs-2));
    sigma  = exp(-tau^2+tau*randn(T+burnin-p,1)); % Second moment = 1
    shocks = zeta.*sigma; % SVAR shocks
    Y      = var_sim(A_true, zeros(n,1), shocks*H_true', T+burnin, zeros(p,n)); % Simulate VAR
    Y      = Y(end-T+1:end,:); % Discard burn-in
    
    % Calculate test statistic and bootstrap quantiles
    the_H_estim = nan(n,n,2);
    [teststats(ii), teststats_boot, ~, ~, ~, the_H_estim(:,:,1), the_H_estim(:,:,2)] ...
        = var_test_indep(Y, p_estim, pml_settings, numboot, false);
    teststats_boot_quants(ii,:) = quantile(teststats_boot,quants);
    
    % Flip signs and permute estimated H to find best fit to truth
    for im=1:2 % For both PML and cumulant estimate...
        the_estim = the_H_estim(:,:,im).*sign(diag(the_H_estim(:,:,im))'); % Flip signs to ensure positive diagonal
        H_estims(ii,:,:,im) = permute_mat(the_estim, H_true); % Best-fitting column permutation
    end
    
    % Print progress
    if mod(ii,ceil(numsim/100))==0
        fprintf('%s%3d%s\n', repmat(' ',1,round(50*ii/numsim)), round(100*ii/numsim), '%');
    end

end

elapsed_time = toc(timer);
disp('Done. Elapsed time (min):');
disp(elapsed_time/60);


%% Display results

disp('True H');
disp(H_true);

disp('Mean estimated H: PML');
disp(squeeze(mean(H_estims(:,:,:,1),1)));

disp('Mean estimated H: cumulant');
disp(squeeze(mean(H_estims(:,:,:,2),1)));

disp('Std. of estimated H: PML');
disp(squeeze(std(H_estims(:,:,:,1),0,1)));

disp('Std. of estimated H: cumulant');
disp(squeeze(std(H_estims(:,:,:,2),0,1)));

disp('Test rejection rate (top: nominal, bottom: actual):');
disp(1-quants);
disp(mean(teststats>teststats_boot_quants));

delete(poolobj);
save(['out_' datestr(now, 'yyyy_mmdd_HHMM') '.mat'], '-v7.3')
