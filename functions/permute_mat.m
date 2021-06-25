function best_fit = permute_mat(mat, target)

% Best-fitting column permutation of matrix "mat" to "target"

column_perms = perms(1:size(mat,2)); % All permutations of column indices
best_fit = mat;
err = Inf;

for ip=1:size(column_perms,1) % Loop over all permutations
    mat_new = mat(:,column_perms(ip,:)); % Permuted matrix
    err_new = norm(mat_new-target, 'fro'); % Error relative to target
    if err_new<err
        best_fit = mat_new;
        err = err_new;
    end
end

end