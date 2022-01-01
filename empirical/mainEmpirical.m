%% Preliminaries
clear
clc
close all

rng(202105272, 'twister');      % Seed RNG


% Settings
runSpec     = 1;      % 0 to load existing results files. 1 to rerun.
vis         = 'off';  % Show figures
pmax        = 8;      % Maximum VAR lag length for lag selection exercise
verbose     = true;
initSetting = 'GlobalSearch';
quants      = [0.68 0.75 0.9 0.95];  % Bootstrap quantiles to store
numboot     = 1000;                   % No. of bootstrap iterations



% PML mixture parameters
pml_settings.ps     = repmat(0.5,2,3);                       % Mixture probabilities
pml_settings.mus    = [repmat( 0.1,1,3); repmat(-0.1,1,3)];  % Mixture means
pml_settings.sigmas = [0.5 1.3 0.7; 1.32 0.54 1.22];         % Mixture std dev's
pml_settings.opts   = optimoptions('fmincon', 'Display', 'notify', ...
    'SpecifyObjectiveGradient',  true, ...
    'SpecifyConstraintGradient', true, ...
    'CheckGradients',            false, ...
    'MaxFunctionEvaluations',    5e3);                       % fmincon options


% Estimation windows
windows = [datetime(1959, 10, 1), datetime(2019, 12,31);   % 1959-2019
    datetime(1985, 1 , 1), datetime(2019, 12 ,31);         % 1985-2019
    datetime(1973, 1,  1), datetime(2019, 12, 31)          % 1973-2019
    ];


df = readtimetable('data/data.xls');  % Read data


% Set paths
pathFunctions   = '../functions/';
pathFiguresBase = 'figures2019/';
pathFigs = [pathFiguresBase 'initSetting=' initSetting '_'];
mkdir(pathFiguresBase)
addpath(pathFunctions)


%% Store empirical specifications of interest

Spec = struct;

% Add estimation windows to the specification file
for jWindow = 1:size(windows, 1)
        
    % Variables and lags
    Spec(jWindow).name      =  'ygap, no controls';
    Spec(jWindow).mnems_VAR = {'pi', 'ygap', 'R'};
    Spec(jWindow).mnems_exo = [];        % Exogenous variables
    Spec(jWindow).p          = 6;         % lags    
        
    % Reference matrix (from GMR 2017)
    Spec(jWindow).CRef = [0.9560   -0.2710   -0.1160;
                          0.2580    0.9590   -0.1160;
                          0.1430    0.0810    0.9860];
    
    
    % Add dates
    Spec(jWindow).TStart = windows(jWindow, 1);
    Spec(jWindow).TEnd   = windows(jWindow, 2);
    
    % Add estimation sample to the specification's name
    Spec(jWindow).nameDate = ...
        [Spec(jWindow).name, ...
        [' (' datestr(windows(jWindow, 1), 'yyyy:qq') '-'  ...
        datestr(windows(jWindow, 2), 'yyyy:qq') ')']];
    
end

nSpec = length(Spec);


%% Run specifications
if runSpec == 0
    load([pathFigs 'Results.mat'])
else
    for j = 1:nSpec
        
        % Store specification
        mnems_VAR = Spec(j).mnems_VAR;
        mnems_exo = Spec(j).mnems_exo;
        p         = Spec(j).p;
        CRef      = Spec(j).CRef;
        TStart    = Spec(j).TStart;
        TEnd      = Spec(j).TEnd;
        dfEst     = df(df.Time >= TStart - calquarters(p) & df.Time <= TEnd, :);
        Y         = dfEst{:, mnems_VAR};  % Endogenous variables
        Z         = dfEst{:, mnems_exo};  % Exogenous variables
        n         = size(Y, 2);
        
        
        % Indexing on margin 3: 1=PML, 2=Cumulant
        C_estim_raw = nan(n, n, 2);  % C from routines
        C_estim     = nan(n, n, 2);  % Permuting to match matrices from paper
        H_estim     = nan(n, n, 2);
        C_eye       = nan(n, n, 2);
        
        % Store lag length objects
        aic         = nan(pmax, 1);
        hq          = nan(pmax, 1);  % Hannan-Quinn
        sc          = nan(pmax, 1);  % Schwarz criterion
        
        
        % Find optimal lag length
        for pTest = 1:pmax
            % Offset by lag length so that sample length is consistent
            dfEst_p     = df(df.Time >= TStart - calquarters(pTest) & df.Time <= TEnd, :);
            Y_p         = dfEst_p{:, mnems_VAR};  % Endogenous variables
            Z_p         = dfEst_p{:, mnems_exo};  % Exogenous variables
            
            [~, Sigma, ~, res, ~, aic(pTest), sc(pTest), hq(pTest)] = ...
                var_estim2(Y_p, pTest, Z_p);            
        end
                
        
        % Run test
        if strcmp(initSetting, 'identity')
            [teststats, teststats_boot, A, c, shocks, ...
                H_estim(:,:,1),     H_estim(:,:,2),...
                C_estim_raw(:,:,1), C_estim_raw(:,:,2), allmins, loglik] = ...
                var_test_indep_emp(Y, Z, p, pml_settings, numboot,...
                                    'verbose', verbose, 'initSetting', eye(n),...
                                    'initSettingBoot', 'full');
            
        else 
            [teststats, teststats_boot, A, c, shocks, ...
                H_estim(:,:,1),         H_estim(:,:,2),...
                C_estim_raw(:,:,1), C_estim_raw(:,:,2), allmins, loglik] = ...
                var_test_indep_emp(Y, Z, p, pml_settings, numboot, ...
                                    'verbose', verbose, 'initSetting', initSetting,...
                                    'initSettingBoot', 'full');                                    
        end
        
        
        teststats_boot_quants = quantile(teststats_boot,quants); 
        
        % Flip signs and permute estimated C to find best fit to reference
        % matrices
        for im=1:2 % For both PML and cumulant estimate
            % Flip signs to ensure positive diagonal
            the_estim       = C_estim_raw(:,:,im).*sign(diag(C_estim_raw(:,:,im))');
            C_estim(:,:,im) = permute_mat(the_estim, CRef); % Best-fitting column permutation
            C_estim_raw(:,:,im) = the_estim;
            
        end
        
        % Store output
        Spec(j).C_estim               = C_estim;
        Spec(j).C_estim_raw           = C_estim_raw;
        Spec(j).H_estim               = H_estim;
        Spec(j).A                     = A;
        Spec(j).c                     = c;
        Spec(j).dfEst                 = dfEst;
        Spec(j).aic                   = aic;
        Spec(j).hq                    = hq;
        Spec(j).sc                    = sc;
        Spec(j).allmins               = allmins;
        Spec(j).teststats             = teststats;
        Spec(j).teststats_boot        = teststats_boot;
        Spec(j).teststats_boot_quants = teststats_boot_quants;
        Spec(j).loglik                = loglik;
        Spec(j).shocks                = shocks;
    end
end

%% Plot data

plotVars = {'pi', 'ugap', 'ygap', 'R'};
f               = figure('Visible', vis);
f.Units         = 'inches';
f.Position(3:4) = [6,4];

for jVar = 1:length(plotVars)
    subplot(3, 2, jVar)
    plot(df.Time, df{:, plotVars{jVar}})
    title(plotVars{jVar})
end

saveas(gcf, [pathFigs 'data.png'])


%% Compare lag lengths for the VAR

idxPlot = 1:3;
Cols = [228,26,28;  % Colors
    55,126,184;
    77,175,74]./255;

f               = figure('Visible', vis);
f.Units         = 'inches';
f.Position(3:4) = [6,6];

plotPosition = 1;

for jFig = idxPlot
    subplot(3,1,plotPosition)
    hold on
    
    % Criteria by lag
    plot(1:pmax, Spec(jFig).aic, 'Color', Cols(1, :))
    plot(1:pmax, Spec(jFig).hq , 'Color', Cols(2, :))
    plot(1:pmax, Spec(jFig).sc , 'Color', Cols(3, :))
    
    % Compute optimal lag length
    [aic_opt, p_aic]  = min(Spec(jFig).aic);
    [hq_opt,  p_hq]   = min(Spec(jFig).hq);
    [sc_opt,  p_sc]   = min(Spec(jFig).sc);
    
    
    % Plot optimal lag length
    scatter(p_aic, aic_opt, 'd', 'filled', ...
        'MarkerEdgeColor', Cols(1, :), 'MarkerFaceColor', Cols(1,:))
    scatter(p_hq, hq_opt, 'd', 'filled', ...
        'MarkerEdgeColor', Cols(2, :), 'MarkerFaceColor', Cols(2,:))
    scatter(p_sc, sc_opt, 'd', 'filled', ...
        'MarkerEdgeColor', Cols(3, :), 'MarkerFaceColor', Cols(3,:))
    title(Spec(jFig).nameDate)
    box on; grid on;
    xlabel('lags')
    legend({'AIC', 'Hannan-Quinn', 'Schwarz'}, 'Location', 'northwest')
    xlim([1 pmax])
    
    plotPosition = plotPosition+1;
end
saveas(gcf, [pathFigs 'optimalLags.png'])



%% Create test rejection rate table

vars_quant = strcat('boot_p', string(100*quants));
out = array2table(nan(nSpec, length(quants)), ...
    'VariableNames', vars_quant);

% Format table
out.specification = string({Spec.nameDate}');
out.teststats     = [Spec.teststats]';
out.pvalue        = nan(nSpec, 1);
out               = out(:, [end-2, end-1, end, 1:end-3]);

for jSpec = 1:nSpec
    teststats              = Spec(jSpec).teststats;
    out{jSpec, vars_quant} = Spec(jSpec).teststats_boot_quants;
    out.pvalue(jSpec)      = mean(Spec(jSpec).teststats_boot > Spec(jSpec).teststats);
end

writetable(out, [pathFigs 'bootTest.xls'])

close all
save([pathFigs 'Results.mat'], '-v7.3')

