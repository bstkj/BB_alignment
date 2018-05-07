function [centroid, n] = best_plane(x, y, z)
centroid = [mean(x); mean(y); mean(z)];
x_t = transpose(x);
y_t = transpose(y);
z_t = transpose(z);
X = vertcat(x_t, y_t, z_t) - centroid;
[U, S, ~] = svd(X);
[~, I] = min(max(S, [], 2));
n = U(:, I);
end

