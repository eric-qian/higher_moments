clear all;
addpath('functions');

% Test ICA code
% MPM 2021-05-26


%% Settings

% DGP
rng(202105261, 'twister');
dfs = [9 11 13];    % df of t-distributions
n = length(dfs);    % No. of series
H_true = randn(n);  % True mixing matrix
T = 1e4;            % Sample size

% PML mixture parameters
pml_settings.ps = repmat(0.5,2,n);                       % Mixture probabilities
pml_settings.mus = [repmat(0.1,1,n); repmat(-0.1,1,n)];  % Mixture means
pml_settings.sigmas = [0.5 1.3 0.7; 1.32 0.54 1.22];     % Mixture std dev's
pml_settings.opts = optimoptions('fmincon', 'Display', 'iter-detailed', ...
                        'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                        'CheckGradients', true); % fmincon options


%% Simulate data

cumuls_true = 6./(dfs-4); % True 4th cumulants of shocks
S = trnd(repmat(dfs,T,1))./sqrt(dfs./(dfs-2)); % Shocks
X = S*H_true'; % Data


%% Estimate

% Estimation by fourth order cumulants
[H_cumul, C_cumul, cumuls, cumul_mat] = ica_cumul(X);

% Estimation by PML
[H_pml, C_pml] = pml(X, @(X) normal_mixture(X, pml_settings.mus, pml_settings.sigmas, pml_settings.ps), C_cumul, pml_settings.opts);
