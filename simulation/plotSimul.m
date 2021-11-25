%% Preliminaries
clear
clc
close all

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

f = figure;
f.Units = 'inches';
f.Position(3:4) = [6,3];

plot(powerTable.df_vol, powerTable.pi_05)
close all



