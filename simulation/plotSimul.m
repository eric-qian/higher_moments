%% Preliminaries
clear
clc
close all

addpath('../functions/')
outDir = 'output/';

load('Results.mat')

%% Organize results

powerTable = array2table(nan(nSpec, 4), 'VariableNames', {'df_vol', 'DGP', 'pi_10', 'pi_05'});
powerTable.DGP = string(powerTable.DGP);

for j = 1:nSpec
    
    powerTable.df_vol(j) = Spec(j).obj.df_vol;
    powerTable{j, 'DGP'}   = string(strrep(func2str(Spec(j).obj.shockDGPSettings.FUN), 'simS_', ''));
    powerTable.pi_10(j) = mean(Spec(j).Sim.pVals < .1);
    powerTable.pi_05(j) = mean(Spec(j).Sim.pVals < .05);
    
end

%% Plot results

plotTableDGP1 = powerTable(strcmp(powerTable.DGP, 'DGP1'), :);
plotTableDGP2 = powerTable(strcmp(powerTable.DGP, 'DGP2'), :);
plotTableDGP3 = powerTable(strcmp(powerTable.DGP, 'DGP3'), :);


%% Make plots

% Plot settings
close all
f = figure;
f.Units = 'inches';
f.Position(3:4) = [6,2];
f.PaperPosition = [0 0 6 2];
f.PaperSize     = [6,2];
hold on

subplot(1,3,1)
plotPower(plotTableDGP1)
ttl = title('Panel A. DGP1');
formatPlot()

subplot(1,3,2)
plotPower(plotTableDGP2)
title('Panel B. DGP2')
formatPlot()

subplot(1,3,3)
plotPower(plotTableDGP3)
title('Panel C. DGP3')
formatPlot()

saveas(gcf, [outDir 'power.pdf'])



%% Compute estimated IRFs


hmax     = 8;
IRF_est  = nan(nSpec, nSim, hmax+1, nVar, nVar);  % sim x hz x var x shock

% Compute IRFs for each specification, simulation, shock
for jSpec = 1:nSpec
    for jSim = 1:nSim
        A_estim = squeeze(Spec(jSpec).Sim.A_estim(jSim, :, :));
        c_estim = squeeze(Spec(jSpec).Sim.c_estim(jSim, :, :));
        H_estim = Spec(jSpec).Sim.H_estim(jSim, :, :);
        [A_c, c_c, H_c] = toCompanion(A_estim, c_estim, H_estim);
        
        for jShock = 1:nVar
            
            % Shock index
            ind         = zeros(nVar, 1);
            ind(jShock) = 1;
            
            IRF_est(jSpec, jSim, :, :, jShock) = getIRF(A_c, H_c, ind, hmax);  % Get IRF
        end
        
    end
end

%% Compute true IRF

IRF_true = nan(hmax+1, nVar, nVar);
[A_true_c, c_true_c, H_true_c] = ...
    toCompanion(obj.VARDGPSettings.A_true, ...
                obj.VARDGPSettings.c_true, ...
                obj.VARDGPSettings.H_true);

            
for jShock = 1:nVar
    
    % Shock index
    ind         = zeros(nVar, 1);
    ind(jShock) = 1;
    
    IRF_true(:, :, jShock) = getIRF(A_true_c, H_true_c, ind, hmax);  % Get IRF
end


%% Plot IRFs

Cols = [228,26,28
    55,126,184
    77,175,74] ./255;

LW = 1.5;
yl = [0 0.5];

size(IRF_est)
size(IRF_true)

IRF_true_mat = permute(repmat(IRF_true, 1, 1, 1, nSpec, nSim), [4, 5, 1:3]);

% Decompose MSE into bias and variance
MSE  = squeeze(mean((IRF_est - IRF_true_mat).^2, 2));
Var  = squeeze(var(IRF_est, 1, 2)) ;
Bias = squeeze(mean(IRF_est - IRF_true_mat, 2));

max(MSE(:) - Var(:) - Bias(:).^2)

% Diagnostics: Plot decomposition for each specification
for jSpec = 1:nSpec
    jShock = 1;
    close all
    f = figure;
    f.Units = 'inches';
    f.Position(3:4) = [6,3];
    p1 = plot(MSE(jSpec, :, 1, jShock), 'LineWidth', LW, 'Color', Cols(1,:));
    hold on
    p2 = plot(Var(jSpec, :, 1, jShock), 'LineWidth', LW, 'Color', Cols(2,:));
    p3 = plot(Bias(jSpec, :, 1, jShock).^2, 'LineWidth', LW, 'Color', Cols(3,:));
    ylim(yl)
    title(jSpec)
    legend([p1, p2, p3], {'MSE', 'Variance', 'Bias^2'}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal')
%    pause
end


%% Plot for selected specifications
% Plot settings

indSpec = [1 1+1*3*14 1+57];
jPlot   = 1;
jShock = 2
close all
f               = figure;
f.Units         = 'inches';
f.Position(3:4) = [6,6];
% MSE
for jSpec = 1:length(indSpec)
    subplot(3, 1, 1)
    hold on
    plot(0:hmax, MSE(indSpec(jSpec), :, 1, jShock), ...
        'LineWidth', LW, 'Color', Cols(jSpec,:))
    title('MSE')
    ylim(yl)
    box on
    grid on
end

% Bias
for jSpec = 1:length(indSpec)
    subplot(3, 1, 2)
    hold on
    p_j = plot(0:hmax, Bias(indSpec(jSpec), :, 1, jShock).^2, ...
        'LineWidth', LW, 'Color', Cols(jSpec,:));
    title('Bias^2')
    ylim(yl)
    box on
    grid on
end

p=[];
% Variance
for jSpec = 1:length(indSpec)
    subplot(3, 1, 3)
    hold on
    p_j = plot(0:hmax, Var(indSpec(jSpec), :, 1, jShock), ...
        'LineWidth', LW, 'Color', Cols(jSpec,:));
    title('Variance')
    ylim(yl)
    box on
    grid on
    p = [p p_j];
end

legend(p, strrep({Spec(indSpec).label}, '_vol', ''), ...
    'Orientation', 'horizontal', 'Location', 'southoutside')


%% Functions

function plotPower(plotTable)


LW = .8;
Cols = [228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0] ./255;


scatter(0, 0.05, 25, '+', 'MarkerEdgeColor', 'b')
hold on
scatter(0, 0.10, 25, '+', 'MarkerEdgeColor', 'b')

p1 = plot(1./plotTable.df_vol, plotTable.pi_05, 'LineWidth', LW, 'Color', Cols(1,:), 'LineStyle', '--');
hold on
p2 = plot(1./plotTable.df_vol, plotTable.pi_10, 'LineWidth', LW, 'LineStyle', '-', 'Color', Cols(1,:));
legend([p1, p2], {'5%', '10%'}, 'Location', 'southeast')
xlabel('1/k')
box on
grid on


end

function formatPlot()
ax = gca;
ax.FontName = 'Times';
ax.TitleHorizontalAlignment = 'left';
ax.FontSize = 8;
end
