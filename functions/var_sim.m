function Y_sim = var_sim(A, c, innovs, T, Y_init)

    % Simulate VAR(p) model
    
    % Dimensions
    [n,np] = size(A);
    p = np/n;
    A_resh = reshape(A,n,n,p);
    
    % Initialize
    Y_sim = [Y_init; zeros(T-p,n)];
    
    % Recursion
    for t=p+1:T
        Y_sim(t,:) = c' + innovs(t-p,:);
        for l=1:p
            Y_sim(t,:) = Y_sim(t,:) + Y_sim(t-l,:)*A_resh(:,:,l)';
        end
    end

end