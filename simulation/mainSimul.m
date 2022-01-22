%% Preliminaries
% Trace out power function by departure from the null
clear
clc
close all

poolobj = parpool('local');
addpath(genpath('../functions/'))


% Simulation settings
quants = [0.68 0.75 0.9 0.95];  % Bootstrap quantiles to store
rng(20220113, 'twister');       % Seed RNG


%  Degrees of freedom for shared volatility. Measures departure from null.
df_vol_vec  = [inf 1./(0.1:.1:1)];

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
SpecEmpirical = load('../empirical/figures2019/initSetting=GlobalSearch_Results.mat', 'Spec');
SpecEmpirical = SpecEmpirical.Spec;

jEmpirical = find(year([SpecEmpirical.TStart]) == 1973 & ...
                  year([SpecEmpirical.TEnd  ]) == 2019);

VARDGPSettings.A_true = SpecEmpirical(jEmpirical).A;
VARDGPSettings.c_true = SpecEmpirical(jEmpirical).c;          
VARDGPSettings.H_true = SpecEmpirical(jEmpirical).H_estim(:,:,1);        
      
clear SpecEmpirical jEmpirical

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
