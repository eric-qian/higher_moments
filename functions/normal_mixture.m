function [log_dens, log_dens_grad] = normal_mixture(X, mus, sigmas, probs)

    % Log densities of bivariate normal mixtures, with gradient

    [T,n] = size(X);
    
    log_dens = nan(T,n);
    log_dens_grad = nan(T,n);
    
    for i=1:n % Loop over data elements
        the_norms = normpdf(X(:,i),mus(:,i)',sigmas(:,i)');
        log_dens(:,i) = log(the_norms*probs(:,i)); % Log density
        if nargout>1
            log_dens_grad(:,i) = exp(-log_dens(:,i)).*((-the_norms.*((X(:,i)-mus(:,i)')./(sigmas(:,i)'.^2)))*probs(:,i)); % Gradient
        end
    end

end