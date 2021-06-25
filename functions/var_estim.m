function [A, Sigma, c, res] = var_estim(Y, p)

    % VAR least-squares estimation

    T = size(Y,1);
    X = [lagmatrix(Y, 1:p) ones(T,1)];
    beta = (X(p+1:end,:)\Y(p+1:end,:))';
    A = beta(:,1:end-1);
    c = beta(:,end);
    res = Y(p+1:end,:)-X(p+1:end,:)*beta';
    Sigma = (res'*res)/(size(X,1)-size(X,2));

end