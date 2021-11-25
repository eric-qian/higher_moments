% plotPower.m  Plot power graph
%% Preliminaries
clear
clc

figPath = 'output/';

% Read data
pValCrossCor      = readtable([figPath 'power_pVal_crossCor.csv']);
pValCrossCor.Var1 = [];

pValDN      = readtable([figPath 'power_pVal_DN.csv']);
pValDN.Var1 = [];

% Plot parameters
blue = [103,169,207]./255; 
red  = [239,138,98] ./255;
LW   = 1.5;

tauVec = 0:.025:1;


%% Get power curves
power = array2table(nan(length(tauVec), 4));
power.Properties.VariableNames = {'crossCor_10', 'DN_10', 'crossCor_5', 'DN_5'};


power.crossCor_10 = mean(pValCrossCor{:,:} < 0.1)';
power.crossCor_5  = mean(pValCrossCor{:,:} < 0.05)';
power.DN_10       = mean(pValDN{:,:} < 0.1)';
power.DN_5        = mean(pValDN{:,:} < 0.05)';

%%
close all
f = figure;
hold on
plot(tauVec, power.crossCor_10, 'Color', blue, 'LineWidth', LW)
plot(tauVec, power.crossCor_5,  'Color', blue, 'LineStyle', '--', 'LineWidth', LW)
plot(tauVec, power.DN_10,       'Color', red, 'LineWidth', LW)
plot(tauVec, power.DN_5,        'Color', red, 'LineStyle', '--', 'LineWidth', LW)
box on; grid on;

ylabel('Rejection rate')
xlabel('\gamma')
f.Units         = 'inches';
f.Position(3:4) = [6,3];
legend({'crossCor 10%', 'crossCor 5%', 'DN 10%', 'DN 5%'}, ...
    'Orientation', 'vertical', 'Location', 'southeast')

saveas(gcf, [figPath 'power.png'])
