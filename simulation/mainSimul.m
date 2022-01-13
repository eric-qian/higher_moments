%% Preliminaries
% Trace out power function by departure from the null
clear
clc
close all

poolobj = parpool('local', 1);
addpath(genpath('../functions/'))

% Simulation settings
quants = [0.68 0.75 0.9 0.95];  % Bootstrap quantiles to store
rng(202105272, 'twister');      % Seed RNG


%  Degrees of freedom for shared volatility. Measures departure from null.
df_vol_vec  = [inf 1./logspace(-2.5, 0, 19)];

% Settings (general)
T       = 200;
burnin  = T/2;
nVar    = 3;
nSim    = 1000;
pMat    = [];
numboot = 500;


% DGP1 settings: Independent t 
shockDGP1Settings     = struct;
shockDGP1Settings.FUN = @simS_DGP1;  
shockDGP1Settings.df  = 5;

% DGP2 settings: Independent Laplace
shockDGP2Settings          = struct;
shockDGP2Settings.FUN      = @simS_DGP2;  
shockDGP2Settings.location = 0;
shockDGP2Settings.scale    = sqrt(1/2);

% DGP3 Settings: Mixture of normals
delta                    = [0.8, 1.2, -1];
kappa                    = [0.06, 0.08, 0.2];
lambda                   = [0.52, 0.4, 0.2];
shockDGP3Settings        = struct;
shockDGP3Settings.FUN    = @simS_DGP3;  
shockDGP3Settings.delta  = delta;
shockDGP3Settings.kappa  = kappa;
shockDGP3Settings.lambda = lambda;

% VAR Settings. Taken from ygap (1973-2019) p=6 empirical specification.
% True VAR(p) matrix [A_1,A_2,...,A_p]. 
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
initSetting     = 'GlobalSearch';
initSettingBoot = 'full';

p_est = 6;  % Number of lags for VAR estimation


%% Create simulation objects for each specification

Spec  = struct;
j     = 1;

for df_vol = df_vol_vec
    Spec(j).obj         = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP1Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);
    Spec(j).label       = ['t, df_vol=' num2str(df_vol)]; 
    Spec(j+1).obj       = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP2Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);
    Spec(j+1).label     = ['Laplace, df_vol=' num2str(df_vol)]; 
    Spec(j+2).obj       = higher_moments_simul(T, burnin, nVar, nSim, pMat, df_vol, numboot, p_est, shockDGP3Settings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot);
    Spec(j+2).label     = ['Mixture, df_vol=' num2str(df_vol)]; 
    j                   = j+3;
end

nSpec = length(Spec);

%% Run simulations


for j = 1:nSpec  % For each shock/volatility spec
    obj                   = Spec(j).obj;
    teststats             = nan(nSim, 1);
    pVals                 = nan(nSim, 1);
    teststats_boot_quants = nan(nSim, length(quants));
    H_estim               = nan(nSim, nVar, nVar); % Last dimension: 1 is PML and 2 is cumulant estimate
    H_estim_cumul         = nan(nSim, nVar, nVar); % Last dimension: 1 is PML and 2 is cumulant estimate
    rand_seeds            = randi(2^32-1,nSim, 1);
    A_estim               = nan([nSim, size(obj.VARDGPSettings.A_true)]);
    c_estim               = nan([nSim, size(obj.VARDGPSettings.c_true)]);
    
    
    tic
    parfor ii = 1:nSim
        
        rng(rand_seeds(ii), 'twister'); % Seed RNG
        
        Res                         = runVARTest(obj);
        pVals(ii)                   = mean(Res.teststats_boot > Res.teststat);
        teststats_boot_quants(ii,:) = quantile(Res.teststats_boot, quants);
        teststats(ii)               = Res.teststat;
        H_estim(ii, :, :)           = Res.H;
        H_estim_cumul(ii, :, :)     = Res.H_cumul;
        A_estim(ii, :, :)           = Res.A;
        c_estim(ii, :, :)           = Res.c;
    
    end
    
    runtime = toc;
    
     % Store simulation objects
     Sim                       = struct;
     Sim.teststats             = teststats;
     Sim.pVals                 = pVals;
     Sim.teststats_boot_quants = teststats_boot_quants;
     Sim.H_estim               = H_estim;
     Sim.H_estim_cumul         = H_estim_cumul;     
     Sim.A_estim               = A_estim;
     Sim.c_estim               = c_estim;
     Sim.rand_seeds            = rand_seeds;
     Sim.runtime               = runtime;
     Spec(j).Sim               = Sim; 
     save('Results.mat', '-v7.3') 
end

delete(poolobj)
