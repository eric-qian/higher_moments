%% Preliminaries
clear
clc
close all

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
close all

plotPower(plotTableDGP1)
title(['DGP1 (T= ' num2str(T) ')'])
saveas(gcf, [outDir 'DGP1.png'])


plotPower(plotTableDGP2)
title(['DGP2 (T= ' num2str(T) ')'])
saveas(gcf, [outDir 'DGP2.png'])

plotPower(plotTableDGP3)
title(['DGP3 (T= ' num2str(T) ')'])
saveas(gcf, [outDir 'DGP3.png'])



function plotPower(plotTable)


LW = 1.25;
Cols = [228,26,28;
    55,126,184;
    77,175,74;
    152,78,163;
    255,127,0] ./255;

f = figure;
f.Units = 'inches';
f.Position(3:4) = [6,3];
hold on
plot(1./plotTable.df_vol, plotTable.pi_05, 'LineWidth', LW, 'Color', Cols(1,:), 'LineStyle', '--')
hold on
plot(1./plotTable.df_vol, plotTable.pi_10, 'LineWidth', LW, 'LineStyle', '-', 'Color', Cols(1,:))
legend({'5%', '10%'}, 'Location', 'southeast')
xlabel('1/df')
box on
grid on
end
