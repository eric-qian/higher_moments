function cumul_mat = cumul4(X)

    % Estimate fourth order cumulant tensor matrix

    [T,n] = size(X); % Data dimensions
    
    % Indices of matrix vectorization
    inds_vec = [repmat((1:n)',n,1) kron((1:n)',ones(n,1))]; % Indices of vec(matrix)
    inds_vech = inds_vec(tril(ones(n))==1,:); % Indices of vech(matrix)
    
    % Expanded data vectors
    Xt1 = X(:,inds_vech(:,1));
    Xt2 = X(:,inds_vech(:,2));
    
    % Cumulant tensor matrix
    aux = Xt1'*Xt2/T;
    cumul_mat = cov(Xt1.*Xt2,1)-(Xt1'*Xt1/T).*(Xt2'*Xt2/T)-aux.*aux';

end