% run_sum_var_runtime.m  Profile code. Written so that it can be run
% locally.

clear
clc
close all
addpath('../functions');
outPath = 'Res/';
mkdir(outPath);

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

% Set tau, numsim, numboot ------------------------------------------------
% Large # of sims/bootstrap draws
TVec    = [50, 100, 250, 500, 1000, 2000, 5000, 10000];                               % Sample size
%TVec    = 2000000;
tau     = 0;
numsim  = 500;  % Could be vector input
numboot = 50;  % No. of bootstrap iterations

% Small # of sims/boot
% tauVec = 0:.01:1;
% numsim  = 1000;
% numboot = 500;

nSpec             = length(TVec);
Spec              = struct;

for jSpec = 1:nSpec
    Spec(jSpec).T             = TVec(jSpec);
end




% Set T and burnin---------------------------------------------------------

% Small sample size length
% T      = 250;                               % Sample size
% burnin = 100;                               % Simulation burn-in


% Set other settings ------------------------------------------------------

% DGP
dfs    = [9 11 13];                         % df of zeta t-distributions
A_true = [0.8 0.3 0; 0 0.7 0; 0.2 0.2 0.2]; % True VAR(p) matrix [A_1,A_2,...,A_p]
H_true = [1 0 0; ones(1,3); 3 2 1];         % True impact impulse responses (diagonal elements should be positive)

% Estimation/test
p_estim = 1;    % Estimation lag length

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
rng(202105272, 'twister');      % Seed RNG

%% % Run specifications

for jSpec = 1:nSpec 
    %% Extract specification-varying objects
    
    T     = Spec(jSpec).T;
    % Long sample length
    burnin = T/2;                               % Simulation burn-in

    
    %% DGP    
    [n,np] = size(A_true);
    p      = np/n;
    
    comp_form = [A_true; eye(np-n) zeros(np-n,n)]; % Companion form
    comp_eigs = eig(comp_form);
    
    disp('Eigenvalues of companion matrix (magnitudes)');
    disp(sort(abs(comp_eigs),'descend'));
    
    
    %% Simulate
    
    teststats             = nan(numsim,1);
    teststats_boot_quants = nan(numsim,length(quants));
    H_estims              = nan(numsim,n,n,2); % Last dimension: 1 is PML and 2 is cumulant estimate
    rand_seeds            = randi(2^32-1,numsim,1);
    
    disp('Simulating...');
    
    elapsed_time       = array2table((1:numsim)', 'VariableNames', {'run'});
    elapsed_time.cumul = nan(numsim, 1);
    elapsed_time.pml   = nan(numsim, 1);    
    
    for ii=1:numsim

        clc
        disp([num2str(ii), ' of ' num2str(numsim) , ' (T=' num2str(T) ')...'])
        rng(rand_seeds(ii), 'twister'); % Seed RNG
        
        % Simulate VAR data
        zeta   = trnd(repmat(dfs,T+burnin-p,1))./sqrt(dfs./(dfs-2));
        sigma  = exp(-tau^2+tau*randn(T+burnin-p,1)); % Second moment = 1
        shocks = zeta.*sigma; % SVAR shocks
        Y      = var_sim(A_true, zeros(n,1), shocks*H_true', T+burnin, zeros(p,n)); % Simulate VAR
        Y      = Y(end-T+1:end,:); % Discard burn-in
        
        % Calculate test statistic and bootstrap quantiles
        the_H_estim = nan(n,n,2);
        [A, ~, c, res] = var_estim(Y, p);
    
        % Preliminary estimate of H from fourth order cumulants
        tic;
        [H_cumul,C_cumul] = ica_cumul(res);
        elapsed_time.cumul(ii) = toc;
        
        tic;
        % Estimate H by PML
        H = pml(res, @(X) normal_mixture(X, pml_settings.mus, pml_settings.sigmas, pml_settings.ps), C_cumul, pml_settings.opts);
        elapsed_time.pml(ii) = toc;        
        
%         [teststats(ii), teststats_boot, ~, ~, ~, the_H_estim(:,:,1), the_H_estim(:,:,2)] ...
%             = var_test_indep(Y, p_estim, pml_settings, numboot, false);
%         teststats_boot_quants(ii,:) = quantile(teststats_boot,quants);
        
%         % Flip signs and permute estimated H to find best fit to truth
%         for im=1:2 % For both PML and cumulant estimate...
%             the_estim = the_H_estim(:,:,im).*sign(diag(the_H_estim(:,:,im))'); % Flip signs to ensure positive diagonal
%             H_estims(ii,:,:,im) = permute_mat(the_estim, H_true); % Best-fitting column permutation
%         end
%         
%         % Print progress
%         if mod(ii,ceil(numsim/100))==0
%             fprintf('%s%3d%s\n', repmat(' ',1,round(50*ii/numsim)), round(100*ii/numsim), '%');
%         end
        
    end
    
    Spec(jSpec).elapsed_time = elapsed_time;
    
    
             
    
end

%% Check runtimes

cumulTime = nan(numsim, nSpec);
pmlTime   = nan(numsim, nSpec);

for jSpec = 1:nSpec
    cumulTime(:, jSpec) = Spec(jSpec).elapsed_time.cumul;
    pmlTime(:, jSpec)   = Spec(jSpec).elapsed_time.pml;    
end

% Densities of cumul times 
close all
figure
subplot(2,1,1)
plot(TVec, sum(pmlTime))
ylabel('seconds')
xlabel('T')
title(['PML total runtime (numsim=' num2str(numsim) ')' ])

subplot(2,1,2)
plot(TVec, sum(cumulTime))
ylabel('seconds')
xlabel('T')
title(['cumul total runtime (numsim=' num2str(numsim) ')' ])
saveas(gcf, 'runtime.png')


