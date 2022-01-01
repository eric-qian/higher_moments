% higher_moments_simul_test.m Tests for the higher_moments_simul
%                             class

%% Preliminaries
clear
clc
close all

addpath(genpath('../functions/'))


TShocks = 1000000;  % Sample size of shocks. Verify result is mean 0 variance 1. 

% Settings (general)
T       = 200;
burnin  = T/2;
nVar    = 3;
nSim    = 100;
pMat    = [];
df_vol  = inf;
numboot = 50;


% DGP1 settings: Independent t 
shockDGP1Settings    = struct;
shockDGP1Settings.FUN = @simS_DGP1;  
shockDGP1Settings.df = 5;

% DGP2 settings: Independent Laplace
shockDGP2Settings          = struct;
shockDGP2Settings.FUN       = @simS_DGP2;  
shockDGP2Settings.location = 0;
shockDGP2Settings.scale    = sqrt(1/2);


% DGP3 Settings: Mixture of normals
delta               = [0.8, 1.2, -1];
kappa               = [0.06, 0.08, 0.2];
lambda              = [0.52, 0.4, 0.2];
shockDGP3Settings        = struct;
shockDGP3Settings.FUN    = @simS_DGP3;  
shockDGP3Settings.delta  = delta;
shockDGP3Settings.kappa  = kappa;
shockDGP3Settings.lambda = lambda;

% VAR Settings
% True VAR(p) matrix [A_1,A_2,...,A_p]. Taken from ygap (1973-2019) p=6
% empirical specification.
VARDGPSettings.A_true = [0.608738397	-0.000866205	 0.014683166;
                      13.81621469	 1.205064592	 33.03418321;
                      0.353920451	 0.002494987	 1.171766543;
                      0.149046271	 0.001737221	 0.213937433;
                     -23.36489153	-0.153184317	-14.88692223;
                     -0.460032721	-0.005666657	-0.74314992;
                     -0.000920109	-0.001652481	-0.066513462;
                      20.19423764	-0.172304161	-20.38470652;
                      0.198278699	 0.004338867	 0.843998171;
                      0.231465008	 0.000464581	-0.032134486;
                      12.76470394	 0.136737418	 3.523379531;
                     -0.160573781	-0.001868048	-0.477656288;
                     -0.104138921	 0.000621191	-0.003086554;
                     -14.54675365	-0.136325608	 6.106542731;
                       0.16933	     0.000215091	 0.38732398;
                       0.069589466	-0.000460019	-0.041969878;
                      -2.553050349	 0.023230699	-6.29112383;
                      -0.115254249	 0.000323461	-0.234394405]';
VARDGPSettings.c_true = [0.278583924, -0.000118978, -0.035718651]';          
 % True impact impulse responses (diagonal elements should be positive)
VARDGPSettings.H_true = [0.267309923	-0.211668089	-0.821604431;
                      0.004657489	-0.002954198	 0.002411566;
                      0.460142403	 0.360429829	-0.030133683];        
      

      
% PML mixture parameters
PMLSettings.ps       = repmat(0.5,2,3);                       % Mixture probabilities
PMLSettings.mus      = [repmat( 0.1,1,3);
                         repmat(-0.1,1,3)];  % Mixture means
PMLSettings.sigmas   = [0.5 1.3 0.7;
    1.32 0.54 1.22];     % Mixture std dev's
PMLSettings.opts     = optimoptions('fmincon', 'Display', 'notify', ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
    'CheckGradients', false, 'MaxFunctionEvaluations', 5e3); % fmincon options      
initSetting = 'GlobalSearch';
initSettingBoot = 'full';

p_est = 6;  % Number of lags for VAR estimation

%% Create simulation objects
obj_t    = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP1Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);
obj_lapl = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP2Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);
obj_mix  = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP3Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);


%% Test t shocks

S   = simS_DGP1(obj_t, TShocks);

% First two moments
cov(S)
mean(S)

%% Test Laplace shocks

S   = simS_DGP2(obj_lapl, TShocks);

% First four moments
cov(S)
mean(S)
skewness(S)
kurtosis(S)

%% Test mixture shocks

S   = simS_DGP3(obj_mix, TShocks);

% First two moments
cov(S)
mean(S)


ksdensity(S(:, 1))
ksdensity(S(:, 2))
ksdensity(S(:, 3))
%% Test stochastic volatility function

S   = simS_vol(obj_t, TShocks);


%% Test data generation function for DGP 1

obj_t.T = 200000;
Y   = simY(obj_t);
[A, Sigma, c] = var_estim2(Y, obj_t.p_est, []);
disp('MSE, A')
disp(mean((A(:)- obj_t.VARDGPSettings.A_true(:)).^2 ))

disp('MSE, c')
disp(mean((c(:)- obj_t.VARDGPSettings.c_true(:)).^2 ))


%% Test data generation function for DGP 3

obj_mix.T = 200000;
Y   = simY(obj_t);
[A, Sigma, c] = var_estim2(Y, obj_mix.p_est, []);
disp('MSE, A')
disp(mean((A(:)- obj_mix.VARDGPSettings.A_true(:)).^2 ))

disp('MSE, c')
disp(mean((c(:)- obj_mix.VARDGPSettings.c_true(:)).^2 ))


%% Test data generation function for short sample, t shocks
obj_t.T = 20000;
Y   = simY(obj_t);
[A, Sigma, c] = var_estim2(Y, obj_t.p_est, []);
disp('MSE, A')
disp(mean((A(:)- obj_t.VARDGPSettings.A_true(:)).^2 ))

disp('MSE, c')
disp(mean((c(:)- obj_t.VARDGPSettings.c_true(:)).^2 ))

%% Test data generation function for short sample, laplace shocks
obj_lapl.T = 20000;
Y   = simY(obj_t);
[A, Sigma, c] = var_estim2(Y, obj_lapl.p_est, []);
disp('MSE, A')
disp(mean((A(:)- obj_lapl.VARDGPSettings.A_true(:)).^2 ))

disp('MSE, c')
disp(mean((c(:)- obj_lapl.VARDGPSettings.c_true(:)).^2 ))



%% Test data generation function for short sample, mixture shocks
obj_mix.T = 20000;
Y   = simY(obj_t);
[A, Sigma, c] = var_estim2(Y, obj_mix.p_est, []);
disp('MSE, A')
disp(mean((A(:)- obj_mix.VARDGPSettings.A_true(:)).^2 ))

disp('MSE, c')
disp(mean((c(:)- obj_mix.VARDGPSettings.c_true(:)).^2 ))

%% Test simulation/independence testing function

obj_t.T = 200;
Res     = runVARTest(obj_t);

