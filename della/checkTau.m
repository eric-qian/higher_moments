%% checkTau.m    See how rejection rate depends on tau. Check numsim=200 and 
% numboot = 500.

clear
clc

ResPath = 'Res/';
%%

tauVec = 0:.01:1;  % Tau of interest

% Loop through specification files
for jSpec = 1:length(tauVec)

    % Load data
    Res = load([ResPath, ...
       'out_numboot=500,numsim=1000,tau=' num2str(tauVec(jSpec)) ',T=250.mat']);
   
   % Get quantiles and make rejection rate object
   if jSpec == 1
      quants = Res.quants;
      RejRate = nan(length(tauVec), length(quants));
   end
   
   % Compute rejection rate
    RejRate(jSpec, :) = mean(Res.teststats>Res.teststats_boot_quants);
end

%% Make figure
close all
figure()
plot(tauVec, RejRate)
lgd = legend(string(1-quants), 'Location', 'southeast');
lgd.Title.String = 'Rejection rate (nominal)';
box on
grid on
title('Actual test rejection rate (numboot=500, numsim=1000)')
xlabel('\tau')





