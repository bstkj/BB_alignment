% -------------------------------------------------------------------------
% [Ben] 04/15/18
%
% -------------------------------------------------------------------------

function [startPt, endPt] = assign_rows(x, y, z)
n = length(x);
dist_mtx = zeros(n, n);
nearest_neighbor = zeros(n, 1);
for i = 1:n
    for j = 1:n
        dist_mtx(i, j) = distance_pts([x(i) y(i) z(i)], [x(j) y(j) z(j)]);
%         dist_mtx(i, j) = norm([x(i) y(i) z(i)] - [x(j) y(j) z(j)]);
    end
    [~, idx] = sort(dist_mtx(i, :), 'ascend');
%     nearest_neighbor(i) = idx(2);
    for k = 2:n
        if norm([x(i) y(i)] - [x(idx(k)) y(idx(k))]) < 5
            nearest_neighbor(i) = idx(k);
            break
        end
    end
end
startPt = transpose(1:n);
endPt = nearest_neighbor;
end