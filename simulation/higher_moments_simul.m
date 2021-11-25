classdef higher_moments_simul
    % higher_moments_simul:   Class for running simulation exercises.
    
    properties
        T           = 200       % VAR sample size
        burnin      = 100       % Burn in
        nVar        = 3         % Number of variables
        nSim        = 1000      % Number of samples to be simulated 
        pMat        = []
        df_vol      = inf       % Shared volatility parameter
        numboot     = 1000      % Number of bootstrap draws for independence test
        p_est       = 1
        shockDGPSettings = []        % DGP for shock process
        VARDGPSettings = []        % Settings for VAR estimation
        PMLSettings = []        % Settings for pseudo ML routine
        initSetting = 'cumul'   % Initialization setting for PML routine
        initSettingBoot = 'MLE'  % Initialiation setting for PML routine in bootstrap
    end
    
    methods
        % Initialize object
        function obj = higher_moments_simul(T, burnin, nVar, nSim, pMat, ...
                df_vol, numboot, p_est, shockDGPSettings, VARDGPSettings, PMLSettings, initSetting, initSettingBoot)
            obj.T                = T;
            obj.burnin           = burnin;
            obj.nVar             = nVar;
            obj.nSim             = nSim;
            obj.pMat             = pMat;
            obj.df_vol           = df_vol;
            obj.numboot          = numboot;
            obj.p_est            = p_est;            
            obj.shockDGPSettings = shockDGPSettings;
            obj.VARDGPSettings   = VARDGPSettings;
            obj.PMLSettings      = PMLSettings;
            obj.initSetting      = initSetting;
            obj.initSettingBoot  = initSettingBoot; 
        end
        
        
        % 3 independent, identical t distributions
        function S = simS_DGP1(obj, TShocks)
            df = obj.shockDGPSettings.df;
            S  = trnd(df, [TShocks, obj.nVar])/ sqrt(df/(df-2));  % Normalize so that sd=1
        end
        
        % 3 independent, identical Laplace distributions
        % X = mu - b sgn(U)log(1-2 |U|) ~ Laplace(mu, b)
        function S = simS_DGP2(obj, TShocks)
            mu          = obj.shockDGPSettings.location;
            b           = obj.shockDGPSettings.scale;
                        
            U = rand([TShocks, obj.nVar]) - 0.5;            
            S = mu - b*sign(U).*log(1-2*abs(U));            
            
        end
        
        
        % 3 DLSMN distributions. Chek the parametrization on page 51 of the appendix for
        % Fiorentini and Sentana (2021).
        function S = simS_DGP3(obj, TShocks)
            
            delta  = obj.shockDGPSettings.delta(:);  % Regulates distance between means
            kappa  = obj.shockDGPSettings.kappa(:);  % Ratio of variances
            lambda = obj.shockDGPSettings.lambda(:); % Controls mixture
            
            % Validate consistency of arguments
            assert(obj.nVar == length(delta))
            assert(obj.nVar == length(kappa))
            assert(obj.nVar == length(lambda))
            
            % Store parameters
            mu_1      = delta.*(1-lambda) ./ sqrt(1+lambda.*(1-lambda).*delta.^2);
            mu_2      = -lambda ./ (1-lambda) .* mu_1;
            sigma2_1  = 1/((1+ lambda.*(1-lambda).*delta.^2).*(lambda + (1-lambda).*kappa));
            sigma2_2  = kappa .* sigma2_1;
            
            % Output matrix
            S = nan(TShocks, obj.nVar);
            
            for j = 1:obj.nVar
                
                % Draw mixture components
                y1   = normrnd(mu_1(j), sqrt(sigma2_1(j)), TShocks, 1);
                y2   = normrnd(mu_2(j), sqrt(sigma2_2(j)), TShocks, 1);
                ind1 = rand(TShocks, 1) < lambda(j);
                
                S(ind1, j)  = y1(ind1);
                S(~ind1, j) = y2(~ind1);
            end
            
            
        end
        
        
        % Run shock simulations. Add shared volatility when df_vol < infinity
        function S = simS_vol(obj, TShocks)
            
            % Generate independent shocks
            S = obj.shockDGPSettings.FUN(obj, TShocks);
            
            % Add shared volatility term
            if obj.df_vol < inf
                vol = gamrnd(obj.df_vol/2, 2, TShocks, 1) ./ obj.df_vol;
                
                %        vol2 = chi2rnd(obj.df_vol, TShocks, 1) ./ obj.df_vol;
                S = S .* vol;
            end
            
            
        end
        
        % Simulate data according to VAR in VARDGPSettings
        function Y = simY(obj)
            % Unpack parameters
            A_true = obj.VARDGPSettings.A_true;
            H_true = obj.VARDGPSettings.H_true;
            c_true = obj.VARDGPSettings.c_true;
            
            % Get lags
            [n,np] = size(A_true);
            p      = np/n;

            
            shocks = simS_vol(obj, obj.T+obj.burnin-p);

            Y = var_sim(A_true, c_true, shocks*H_true', obj.T+obj.burnin, zeros(p, obj.nVar)); % Simulate VAR
            Y = Y(end-obj.T+1:end,:); % Discard burn-in
        
        end
        
        % Simulate data and run independence test
        function Res = runVARTest(obj)
            
            Res = struct;
            
            Y = simY(obj);
            [teststat, teststats_boot, A, c, ~, ...
                H, H_cumul, C, C_cumul] = ...
                var_test_indep_emp(Y, [], obj.p_est, obj.PMLSettings, ...
                obj.numboot, 0, obj.initSetting, obj.initSettingBoot);
            
            % Best-fitting permutation
            H_cumul = higher_moments_simul.permuteH(H_cumul, obj.VARDGPSettings.H_true);
            H       = higher_moments_simul.permuteH(H,       obj.VARDGPSettings.H_true);
            
            % Store results
            Res.teststat       = teststat;
            Res.teststats_boot = teststats_boot;
            Res.A              = A;
            Res.c              = c;
            Res.H              = H;
            Res.H_cumul        = H_cumul;
            Res.C              = C;
            Res.C_cumul        = C_cumul;
            
        end
        
    end
    methods (Static)
        function H_perm = permuteH(H, H_true)

            temp   = H .* sign(diag(H)'); % Flip signs to ensure positive diagonal
            H_perm = permute_mat(temp, H_true); % Best-fitting column permutation


        end
    end
end

