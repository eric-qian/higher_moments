function [H, C_opt, allmins, loglik] = pml(X, log_dens, C_init, varargin)
%% Pseudo Maximum Likelihood (PML) estimation of ICA model
% Input:
% - X:             Data.
% - log_dens:      Log density function (first argument is the log density 
%                  and second argument is the gradient).
% - C_init:        Initial whitened mixing matrix.
% - opts:          Optimization settings.
% - globalSetting: 'off' to run fmincon starting at C_init. 'GlobalSearch'
%                  to run the GlobalSearch algorithm.
%
% Output:
% - H:             Optimal unwhitened mixing matrix.
% - C_opt:         Optimal whitened mixing matrix.
% - allmins:       If global optimization routine is run, stores ALL
%                  values of C.
% - loglik:        Log-likelihood.


%% Parse parameters

% Default settings
opts_default          = optimoptions('fmincon', 'Display', 'off', ...
        'SpecifyObjectiveGradient', true, ...
        'SpecifyConstraintGradient', true);
globalSetting_default = 'GlobalSearch';

    
% Setup parser    
parser = inputParser;
addRequired(parser, 'X', @isnumeric);
addRequired(parser, 'log_dens');
addRequired(parser, 'C_init',        @isnumeric);
addOptional(parser, 'opts',          opts_default);
addOptional(parser, 'globalSetting', globalSetting_default, @(x)(ischar(x) | isstring(x)));


parse(parser, X, log_dens, C_init, varargin{:});
opts          = parser.Results.opts;
globalSetting = parser.Results.globalSetting;



% Whiten the data
[T,n] = size(X);
[U,Sigma_chol] = whiten(X);

% Optimization settings
if isempty(C_init)
    C_init = eye(n);
end


% Compute pseudo-MLE of vec(C)
nll = @(C_vec) obj(C_vec,log_dens,U,T,n); % Negative log likelihood

if strcmp(globalSetting, "off")
    [C_vec_opt,  f_opt]= fmincon(nll, C_init(:), ...
        [], [], [], [], [], [], ...
        @(C_vec) cons(C_vec,n), ...
        opts);
    
    allmins = [];
    
elseif strcmp(globalSetting, "GlobalSearch")  % Global search
    gs      = GlobalSearch;
    nParams = length(C_init(:));
    
    problem = createOptimProblem('fmincon', 'objective', nll, ...
        'x0', C_init(:), ...
        'lb', repmat(-1.01, nParams, 1),...  % C is an orthogonal matrix
        'ub', repmat(1.01, nParams, 1), 'options', opts, ...
        'nonlcon', @(C_vec) cons(C_vec,n));
    [C_vec_opt, f_opt,~,~,allmins] = run(gs, problem);
        
    
else
    error('Enter valid globalSetting parameter...')
    
end

loglik = -f_opt;
C_opt  = reshape(C_vec_opt,n,n);  % vec(C) -> C
H      = Sigma_chol*C_opt;        % Un-whiten estimate

end

%% (Negative) log likelihood and gradient

function [val, grad] = obj(C_vec, log_dens, U, T, n)


C = reshape(C_vec,n,n); % vec(C) -> C
CpU = U*C;

if nargout>1
    [the_log_dens,the_log_dens_grad] = log_dens(CpU);
else
    the_log_dens = log_dens(CpU);
end
val = -sum(the_log_dens(:))/T; % Log likelihood

if nargout>1 % If gradient is requested...
    grad = -(U'*the_log_dens_grad)/T; % Gradient wrt. C
    grad = grad(:); % Vectorize
end

end

%% Constraint function 
function [cons_ineq, cons_eq, cons_ineq_grad, cons_eq_grad] = cons(C_vec, n)

% Orthogonality constraint function and gradient

cons_ineq = []; % No inequality constraints

% Equality constraint: C'*C=I
C = reshape(C_vec,n,n); % vec(C) -> C
cons_eq = C'*C-eye(n);
inds_vech = (tril(ones(n))==1); % Only use lower diagonal elements in equation system due to symmetry
cons_eq = cons_eq(inds_vech(:));

if nargout>2 % If gradient is requested...
    % See "Notes on Matrix Calculus" by Paul L. Fackler
    cons_ineq_grad = [];
    aux = kron(eye(n),C');
    inds = reshape(reshape(1:n^2,n,n)',n^2,1); % Indices that transform vec(A) -> vec(A')
    cons_eq_grad = aux(inds_vech,:) + aux(inds(inds_vech(:)),:); % Jacobian
    cons_eq_grad = cons_eq_grad'; % Matlab wants the transpose of the Jacobian
end

end