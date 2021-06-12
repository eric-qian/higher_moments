function [H,C,cumuls,cumul_mat] = ica_cumul(X)

    % ICA using fourth order cumulants

    % Whiten the data
    n = size(X,2);
    [U,Sigma_chol] = whiten(X);
    
    % Estimate cumulant tensor matrix
    cumul_mat = cumul4(U);
    
    % Eigen-decomposition of cumulant tensor matrix
    [V,D] = eig(cumul_mat);
    [eigs_sort,I] = sort(diag(D), 'descend');
    cumuls = eigs_sort(1:n); % n largest eigenvalues
    V_large = V(:,I(1:n));
    
    % Estimate mixing matrix
    C = nan(n);
    low_tri = (tril(ones(n))==1);
    
    for i=1:n
        
        % Construct matrix from vech elements
        the_M = nan(n);
        the_M(low_tri) = V_large(:,i);
        the_M_transpose = the_M';
        the_M(~low_tri) = the_M_transpose(~low_tri);
        
        % Largest eigenvector
        [C(:,i),~] = eigs(the_M,1);
        
    end
    
    % Un-whiten the estimate
    H = Sigma_chol*C;

end