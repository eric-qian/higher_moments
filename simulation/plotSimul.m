%% Preliminaries
clear
clc
close all

addpath('../functions/')
outDir = 'output/';

load('Results.mat')

%% Organize results

powerTable     = array2table(nan(nSpec, 4), ...
                    'VariableNames', {'df_vol', 'DGP', 'pi_10', 'pi_05'});
powerTable.DGP = string(powerTable.DGP);

for j = 1:nSpec
    
    powerTable.df_vol(j) = Spec(j).obj.df_vol;
    powerTable{j, 'DGP'} = string(strrep(func2str(Spec(j).obj.shockDGPSettings.FUN), 'simS_', ''));
    powerTable.pi_10(j)  = mean(Spec(j).Sim.pVals < .1);
    powerTable.pi_05(j)  = mean(Spec(j).Sim.pVals < .05);
    
end

%% Plot results

plotTableDGP1 = powerTable(strcmp(powerTable.DGP, 'DGP1'), :);
plotTableDGP2 = powerTable(strcmp(powerTable.DGP, 'DGP2'), :);
plotTableDGP3 = powerTable(strcmp(powerTable.DGP, 'DGP3'), :);


%% Make plots

% Plot settings
close all
f               = figure;
f.Units         = 'inches';
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




%% Functions

% Formatting for creating power plot
function plotPower(plotTable)

LW = .8;  % Line width
Cols = [228,26,28;  % Colors
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0] ./255;

% Reference markers for 0.05 and 0.10
scatter(0, 0.05, 25, '+', 'MarkerEdgeColor', 'b')
hold on
scatter(0, 0.10, 25, '+', 'MarkerEdgeColor', 'b')

% Plot curves 
p1 = plot(1./plotTable.df_vol, plotTable.pi_05, 'LineWidth', LW, ...
    'Color', Cols(1,:), 'LineStyle', '--');
hold on
p2 = plot(1./plotTable.df_vol, plotTable.pi_10, 'LineWidth', LW, ...
    'LineStyle', '-', 'Color', Cols(1,:));
legend([p1, p2], {'5%', '10%'}, 'Location', 'southeast')

% Add labeling, box, and grid
xlabel('1/k')
box on
grid on


end


% Formatting for plot area
function formatPlot()
ax                          = gca;
ax.FontName                 = 'Times';
ax.TitleHorizontalAlignment = 'left';
ax.FontSize                 = 8;
end
