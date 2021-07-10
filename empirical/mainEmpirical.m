%% Preliminaries
clear
clc
close all

pathFunctions = '../functions/';
addpath(pathFunctions)

runSpec = 0;      % 0 to load existing results files. 1 to rerun.
vis     = 'off';  % Show figures

df = readtimetable('data/data.xls');  % Read data

% Estimation windows
windows = [datetime(1959, 10, 1), datetime(2015, 1, 1);   % Largest sample (1959-2015)
    datetime(1959, 10, 1), datetime(1984, 10,1);          % Pre 84
    datetime(1985, 1 , 1), datetime(2015, 1 ,1);          % Post 84 (1985-2015)
    datetime(1973, 1,  1), datetime(2015, 1, 1)           % Post Bretton Woods (1973 2015)
    ];


hmax        = 48;  % Maximum horizons for IRFs
pmax        = 8;   % Maximum VAR lag length for lag selection exercise

%initSetting = 'cumul';          % Fourth-order cumulant
%initSetting = 'GlobalSearch';  % Use GlobalSearch algorithm
%initSetting = 'MultiStart';  % Use GlobalSearch algorithm
%initSetting = 'identity';      % Identity matrix
%initSetting = 'reference';     % From paper

initSettingVec = {'cumul', 'identity', 'reference', 'GlobalSearch', 'MultiStart'};



for jInit = 1:length(initSettingVec)
    % Initialization settings for PML routine
    
    disp(['Spec ' num2str(jInit) ' of ' num2str(length(initSettingVec)) '...'])
    initSetting = initSettingVec{jInit};
    
    pathFigs = ['figures/initSetting=' initSetting '_'];
    
    
    
    %% Store empirical settings specifications
    % Will run specifications for [variables] x [estimation periods]. Begin by
    % defining variables.
    
    Spec_noest = struct;  % Specification object *without* estimation window
    
    % Unemployment gap no controls
    j = 1;
    Spec_noest(j).name      =  'ugap, no controls';  % Specification name
    Spec_noest(j).mnems_VAR = {'pi', 'ugap', 'R'};   % VAR variables
    Spec_noest(j).mnems_exo = [];                    % Exogenous variables
    
    % Unemployment gap, oil control
    j = 2;
    Spec_noest(j).name      =  'ugap, oil';
    Spec_noest(j).mnems_VAR = {'pi', 'ugap', 'R'};
    Spec_noest(j).mnems_exo = {'doil'};
    
    % Output gap, no controls
    j = 3;
    Spec_noest(j).name      =  'ygap, no controls';
    Spec_noest(j).mnems_VAR = {'pi', 'ygap', 'R'};
    Spec_noest(j).mnems_exo = [];
    
    
    % Output gap, oil control
    j = 4;
    Spec_noest(j).name      =  'ygap, oil';
    Spec_noest(j).mnems_VAR = {'pi', 'ygap', 'R'};
    Spec_noest(j).mnems_exo = {'doil'};
    
    % Other settings
    [Spec_noest.p]        = deal(6);   % lags
    [Spec_noest.TStart]   = deal([]);  % Start date
    [Spec_noest.TEnd]     = deal([]);  % End date
    [Spec_noest.nameDate] = deal([]);  % Name + date
    
    % Unemployment gap reference matrix
    [Spec_noest(1:2).CRef] = deal(    [0.9560   -0.2710   -0.1160;
        0.2580    0.9590   -0.1160;
        0.1430    0.0810    0.9860]);
    
    % Output gap reference matrix
    [Spec_noest(3:4).CRef] = deal([     0.9440    0.3210   -0.0750;
        -0.3270    0.9400   -0.0990;
        0.0390    0.1180    0.9920]);
    
    %% Add estimation windows to the specification object
    
    Spec = Spec_noest;
    
    % Add estimation windows to the specification file
    for jWindow = 1:size(windows, 1)
        
        if jWindow > 1
            Spec = [Spec Spec_noest];
        end
        
        % Index for specifications of interest
        idx = ((jWindow - 1)*4 + 1):((jWindow - 1)*4 + 4);
        
        
        for jj = idx(:)'  % Loop through variables
            
            % Add dates
            Spec(jj).TStart = windows(jWindow, 1);
            Spec(jj).TEnd   = windows(jWindow, 2);
            
            % Add estimation sample to the specification's name
            Spec(jj).nameDate = [Spec(jj).name, ...
                [' (' datestr(windows(jWindow, 1), 'yyyy:qq') '-'  ...
                datestr(windows(jWindow, 2), 'yyyy:qq') ')']];
        end
        
    end
    
    nSpec = length(Spec);
    
    %% Estimation settings
    
    % PML mixture parameters
    pml_settings.ps     = repmat(0.5,2,3);                       % Mixture probabilities
    pml_settings.mus    = [repmat( 0.1,1,3);
        repmat(-0.1,1,3)];  % Mixture means
    pml_settings.sigmas = [0.5 1.3 0.7;
        1.32 0.54 1.22];     % Mixture std dev's
    pml_settings.opts   = optimoptions('fmincon', 'Display', 'notify', ...
        'SpecifyObjectiveGradient', true, ...
        'SpecifyConstraintGradient', true, ...
        'CheckGradients', false, ...
        'MaxFunctionEvaluations', 5e3); % fmincon options
    
    
    % Simulation settings
    quants  = [0.68 0.75 0.9 0.95];  % Bootstrap quantiles to store
    
    % Don't run bootstrap procedure for the GlobalSearch and MultiStart cases
    if strcmp(initSetting, 'GlobalSearch') || strcmp(initSetting, 'MultiStart')
        numboot = [];
    else
        numboot = 1000;                   % No. of bootstrap iterations
    end
    
    verbose = true;
    
    rng(202105272, 'twister');      % Seed RNG
    
    runSpec = 0;
    
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
            
            
            if strcmp(initSetting, 'identity')
                [teststats, teststats_boot, A, c, shocks, ...
                    H_estim(:,:,1),         H_estim(:,:,2),...
                    C_estim_raw(:,:,1), C_estim_raw(:,:,2), allmins, loglik] = ...
                    var_test_indep_emp(Y, Z, p, pml_settings, numboot, verbose, eye(n));
                
            elseif strcmp(initSetting, 'reference')
                [teststats, teststats_boot, A, c, shocks, ...
                    H_estim(:,:,1),         H_estim(:,:,2),...
                    C_estim_raw(:,:,1), C_estim_raw(:,:,2), allmins, loglik] = ...
                    var_test_indep_emp(Y, Z, p, pml_settings, numboot, verbose, CRef);
            else
                [teststats, teststats_boot, A, c, shocks, ...
                    H_estim(:,:,1),         H_estim(:,:,2),...
                    C_estim_raw(:,:,1), C_estim_raw(:,:,2), allmins, loglik] = ...
                    var_test_indep_emp(Y, Z, p, pml_settings, numboot, verbose, initSetting);
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
            Spec(j).C_estim        = C_estim;
            Spec(j).C_estim_raw    = C_estim_raw;
            Spec(j).H_estim        = H_estim;
            Spec(j).dfEst          = dfEst;
            Spec(j).aic            = aic;
            Spec(j).hq             = hq;
            Spec(j).sc             = sc;
            Spec(j).allmins        = allmins;
            Spec(j).teststats      = teststats;
            Spec(j).teststats_boot = teststats_boot;
            Spec(j).teststats_boot_quants = teststats_boot_quants;
            Spec(j).loglik         = loglik;
        end
    end
    
    %% Plot data
    
    plotVars = {'pi', 'ugap', 'ygap', 'R', 'doil'};
    f               = figure('Visible', vis);
    f.Units         = 'inches';
    f.Position(3:4) = [6,4];
    
    for jVar = 1:length(plotVars)
        subplot(3, 2, jVar)
        plot(df.Time, df{:, plotVars{jVar}})
        title(plotVars{jVar})
    end
    
    saveas(gcf, [pathFigs 'data.png'])
    
    %% Output C matrices
    
    % Make red-white-green colormap for easier viewing
    green  = [44,162,95] ./255;
    red    = [222,45,38] ./ 255;
    nColor = 100;
    greens = [linspace(red(1), 1, nColor)', linspace(red(2), 1, nColor)', linspace(red(3), 1, nColor)'];
    reds   = flipud([linspace(green(1), 1, nColor)', linspace(green(2), 1, nColor)', linspace(green(3), 1, nColor)']);
    cm     = [greens; reds];
    
    % Store indices of specifications of interest
    idx_ugap = find(contains({Spec.name}, 'ugap'));
    idx_ygap = find(contains({Spec.name}, 'ygap'));
    
    
    % Compare C matrix for ugap -----------------------------------------------
    figure('Visible', vis)
    f               = gcf;
    f.Units         = 'inches';
    f.Position(3:4) = [10, 6];
    figIdx          = 1;
    
    % Loop through ugap oil figures
    for jFig = idx_ugap
        subplot(3,3, figIdx)
        heatmap(Spec(jFig).C_estim(:, :, 1), 'ColorbarVisible', 'off', 'ColorLimits', [-1 1])
        title(Spec(jFig).nameDate)
        colormap(cm)
        figIdx = figIdx + 1;
    end
    subplot(3,3,9)
    heatmap(Spec(jFig).CRef, 'ColorbarVisible', 'on', 'ColorLimits', [-1 1])
    title('Paper')
    colormap(cm)
    sgtitle(['C, ugap (initSetting=' initSetting ')'])
    saveas(gcf, [pathFigs 'C_ugap_initSetting=' num2str(initSetting) '.png'])
    
    
    % Compare C matrix for ygap -----------------------------------------------
    figure('Visible', vis)
    f = gcf;
    f.Units = 'inches';
    f.Position(3:4) = [10, 6];
    figIdx = 1;
    
    % Loop through ygap oil figures
    for jFig = idx_ygap
        subplot(3,3, figIdx)
        heatmap(Spec(jFig).C_estim(:, :, 1), 'ColorbarVisible', 'off', 'ColorLimits', [-1 1])
        title(Spec(jFig).nameDate)
        colormap(cm)
        figIdx = figIdx + 1;
    end
    subplot(3,3,9)
    heatmap(Spec(jFig).CRef, 'ColorbarVisible', 'on', 'ColorLimits', [-1 1])
    title('Paper')
    colormap(cm)
    sgtitle(['C, ygap (initSetting=' initSetting ')'])
    saveas(gcf, [pathFigs 'C_ygap_initSetting=' num2str(initSetting) '.png'])
    
    
    %% Validation: Compare lag lengths for the VAR
    
    idxPlot = [2, 4];
    Cols = [228,26,28;
        55,126,184;
        77,175,74]./255;
    
    f               = figure('Visible', vis);
    f.Units         = 'inches';
    f.Position(3:4) = [6,6];
    
    plotPosition = 1;
    
    for jFig = idxPlot
        subplot(2,1,plotPosition)
        hold on
        
        % Criteria by lag
        plot(1:pmax, Spec(jFig).aic, 'Color', Cols(1, :))
        plot(1:pmax, Spec(jFig).hq , 'Color', Cols(2, :))
        plot(1:pmax, Spec(jFig).sc , 'Color', Cols(3, :))
        
        % Compute optimal lag length
        [aic_opt, p_aic] = min(Spec(jFig).aic);
        [hq_opt, p_hq]   = min(Spec(jFig).hq);
        [sc_opt, p_sc]   = min(Spec(jFig).sc);
        
        
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
    
    
    %% Check other optima
    
    if strcmp(initSetting, 'GlobalSearch') || strcmp(initSetting, 'MultiStart')
        mkdir([pathFigs, 'checkOptima/'])
        
        for j = 1:nSpec
            allmins        = Spec(j).allmins;
            XOptVals       = [allmins.X];
            FOptVals       = [allmins.Fval];
            [FOpt, jOpt]   = min(FOptVals);
            
            normXOptVals = vecnorm(XOptVals - repmat(XOptVals(:, jOpt), 1, size(XOptVals, 2)));
            
            figure('Visible', vis)
            scatter(normXOptVals, FOptVals, 'filled')
            box on
            grid on
            xlabel('Euclidean distance from optimal value')
            ylabel('-log likelihood')
            saveas(gcf, [pathFigs 'checkOptima/checkOptima_' strrep(Spec(j).nameDate, ':', '') '.png'])
            
            close all
            for jVar1 = 1:size(XOptVals, 1)
                
                close
                f               = figure('Visible', vis);
                f.Units         = 'inches';
                f.Position(3:4) = [6,6];
                figIdx          = 1;
                for jVar2 = 1:size(XOptVals, 1)
                    
                    subplot(3,3, figIdx)
                    [FOptValsSorted, jSorted] = sort(FOptVals, 'descend');
                    p = scatter(XOptVals(jVar1, jSorted), XOptVals(jVar2,jSorted), [], ...
                        FOptValsSorted, 'filled');
                    
                    
                    p.MarkerFaceAlpha = .6;
                    p.SizeData = 20;
                    colormap(gca, 'hot')
                    box on
                    grid on
                    figIdx = figIdx+1;
                    xlabel(['var ' num2str(jVar1)])
                    ylabel(['var ' num2str(jVar2)])
                    hold on
                    scatter(XOptVals(jVar1, jOpt), XOptVals(jVar2, jOpt), ...
                        [], FOpt, '^', 'filled', ...
                        'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
                end
                sgtitle(['Variable ' num2str(jVar1)])
                saveas(gcf, [pathFigs 'checkOptima/plotOptima_' strrep(Spec(j).nameDate, ':', '')...
                    '_var=' num2str(jVar1) '.png'])
                
            end
            
        end
    end
    
    
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
end


%% Summarize results for initialization

loglik_mat = nan(nSpec, length(initSettingVec));
aad_mat    = nan(nSpec, length(initSettingVec));


for jSpec = 1:nSpec
    for jInit = 1:length(initSettingVec)
        % Load results file
        pathFigs = ['figures/initSetting=' initSettingVec{jInit} '_'];
        Results  = load([pathFigs 'Results.mat']);
        Spec     = Results.Spec;
        T        = height(Spec(jSpec).dfEst) - Spec(jSpec).p;  % Sample length
        
        % Store C matrices
        C_pml = Spec(jSpec).C_estim(:,:,1);
        CRef  = Spec(jSpec).CRef; 
        
        loglik_mat(jSpec, jInit) = Spec(jSpec).loglik*T;  % Scaled for comparability

        aad_mat(jSpec, jInit)    = mean(abs(CRef(:) - C_pml(:)));
        
    end
    
    
end


% Show log likelihoods
close all
f = figure('Visible', vis);
h              = heatmap( initSettingVec, {Spec.nameDate}, loglik_mat);
h.ColorScaling = 'scaledrows';
title('Log likelihood')
h.ColorbarVisible = 'off';
f.Units = 'inches';
f.Position(3:4) = [8,8];
saveas(gcf, 'figures/loglik.png')

% Show AAD
close all
f = figure('Visible', vis);
h              = heatmap( initSettingVec, {Spec.nameDate}, aad_mat);
h.ColorScaling = 'scaledrows';
title('Average absolute deviation from GMR')
h.ColorbarVisible = 'off';
f.Units = 'inches';
f.Position(3:4) = [8,8];
saveas(gcf, 'figures/aad.png')


% Output matrices to Excel
out = num2cell(loglik_mat);
out = [{Spec.nameDate}' out];
out = [{''} initSettingVec; out];
writecell(out, 'figures/loglik.xls')

out = num2cell(aad_mat);
out = [{Spec.nameDate}' out];
out = [{''} initSettingVec; out];
writecell(out, 'figures/aad.xls')


