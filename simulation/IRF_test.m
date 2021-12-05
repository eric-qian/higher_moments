%% Preliminaries
% Validate IRF. Example is from Lutkepohl (2005) chapter 3.2.2
clear
clc
close all

hmax = 12;
% 

% 
% c_true = [0.278583924, -0.000118978, -0.035718651]';

A_true = [ -.320 .146 .961 -.161 .115 .934;
            .044 -.153 .289 .050 .019 -.010;
           -.002 .225 -.264 .034 .355 -.022];
c_true = [-.017, .016, .013]';

vcov = [2.130e-03 7.162e-05 1.232e-04;
7.162e-05 1.373e-04 6.146e-05;
1.232e-04 6.146e-05 8.920e-05];
H_true = chol(vcov, 'lower');
H_true = eye(3);


% True impact impulse responses (diagonal elements should be positive)
% H_true = [0.267309923	-0.211668089	-0.821604431;
%     0.004657489	-0.002954198	 0.002411566;
%     0.460142403	 0.360429829	-0.030133683];

n   = size(A_true, 1);
IRF = nan(hmax+1, n, n);  % horizon x variable x shock

%% Get IRFs

[A_c, c_c, H_c] = toCompanion(A_true, c_true, H_true);

for jShock = 1:n

    % Shock index
    ind         = zeros(n, 1);
    ind(jShock) = 1;

    IRF(:, :, jShock) = getIRF(A_c, H_c, ind, hmax);  % Get IRF
end


% Plot IRFs
f = figure();
f.Units = 'inches';
f.Position(3:4) = [6,4];

jPlot = 1;

for jShock = 1:n    
    for jVar = 1:n
        subplot(n, n,jPlot)
        plot(0:hmax, IRF(:, jVar, jShock), 'LineWidth', 1.5)
        hold on
        plot([0, hmax], [0 0], 'Color', 'k', 'LineWidth', 1)
        
        
        title(['Resp of ' num2str(jVar) ' to ' num2str(jShock)])
        jPlot = jPlot + 1;
    end
end

%% Functions
