function [H, C_opt] = pml(X, log_dens, C_init, opts)

    % Pseudo Maximum Likelihood (PML) estimation of ICA model

    % Whiten the data
    [T,n] = size(X);
    [U,Sigma_chol] = whiten(X);
    
    % Optimization settings
    if isempty(C_init)
        C_init = eye(n);
    end
    
    if isempty(opts)
        opts = optimoptions('fmincon', 'Display', 'off', ...
                            'SpecifyObjectiveGradient', true, ...
                            'SpecifyConstraintGradient', true);
    end
    
    % Compute pseudo-MLE of vec(C)
    nll = @(C_vec) obj(C_vec,log_dens,U,T,n); % Negative log likelihood
    C_vec_opt = fmincon(nll, C_init(:), ...
                        [], [], [], [], [], [], ...
                        @(C_vec) cons(C_vec,n), ...
                        opts);
    C_opt = reshape(C_vec_opt,n,n); % vec(C) -> C
    
    % Un-whiten estimate
    H = Sigma_chol*C_opt;

end


function [val, grad] = obj(C_vec, log_dens, U, T, n)

    % (Negative) log likelihood and gradient

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