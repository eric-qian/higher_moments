function [A_c, c_c, H_c] = toCompanion(A_true, c_true, H_true)
p = size(A_true, 2) / size(A_true, 1);
n = size(A_true, 1);

% Companion form
A_c = zeros(p*n);
c_c = zeros(p*n, 1);
H_c = zeros(p*n);

A_c(1:n, :)           = A_true;
A_c(n+1:end, 1:end-n) = eye(n*(p-1));
c_c(1:n)              = c_true;
H_c(1:n, 1:n)         = H_true;
end