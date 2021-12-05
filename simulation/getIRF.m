function IRF = getIRF(A_c, H_c, ind, hmax)

% Get constants
n = length(ind);
p = size(A_c, 2) / n;


shock      = zeros(n*p, 1);
shock(1:n) = ind;

% Add 1 for h = 0
IRF            = nan(hmax+1, n);  
iota           = zeros(n, n*p);
iota(1:n, 1:n) = eye(n); % selection matrix: get non-lagged elements


for h = 0:hmax
    IRF(h+1, :) = iota * A_c^h*H_c*shock;
end

end