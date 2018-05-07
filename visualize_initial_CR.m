% -------------------------------------------------------------------------
% [Ben] 01/16/18
% 3D/2D visualization of the initial ciliary row assignment of basal bodies
% of a cell for given .nd2 file in /nd2_files. This function takes in a
% path to a .nd2 file, and takes about 16s to run. Depending on what the
% parameter 'slice' is set to, either a 3D scatter plot or some 2D slice of
% that 3D plot is produced. slice can be set to 'XYZ', 'XY', 'XZ', or 'YZ'.
% For example: visualize_initial_CR('nd2_files/001.nd2', 'XYZ')
% -------------------------------------------------------------------------

function visualize_initial_CR(imagepath, slice)
channel = 2; % legacy settings
th = 6; % legacy settings
I = readBioImg(imagepath, channel);

[cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, ~, ~, ~, ~] = getBBIdx(I, th);
x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);

OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);
% [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z] = find_cortBB_OAregion(I, th);
% OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);

[pole1, pole2] = findPoles(cort_x, cort_y, cort_z);
d1 = distance_pts(pole1, OA);
d2 = distance_pts(pole2, OA);
% anterior pole is closer to OA region
if d1 < d2
    antPole = pole1;
    postPole = pole2;
else
    antPole = pole2;
    postPole = pole1;
end

num_cort = length(cort_x);
dist2Ant = zeros(num_cort, 1);
for i = 1:num_cort
    dist2Ant(i) = distance_pts([cort_x(i) cort_y(i) cort_z(i)], antPole);
end

[startPt, endPt] = findLink(cort_x, cort_y, cort_z, antPole, postPole, dist2Ant);
% matrix showing ciliary rows
OriginalBBsMatrix = link(startPt, endPt);

% % withLabel is a column vector, withoutLabel is a row vector
% withLabel = intersect(1:num_cort, OriginalBBsMatrix);
% withoutLabel = setdiff(1:num_cort, withLabel);
% minBBsInRow = 10; % legacy settings
% [~, updatedTraceback] = getIniLabel(OriginalBBsMatrix, withLabel,...
%                                     withoutLabel, dist2Ant, ...
%                                     cort_x, cort_y, cort_z, ...
%                                     antPole, postPole, minBBsInRow);
% OriginalBBsMatrix = updatedTraceback;
                               
[num_rows, ~] = size(OriginalBBsMatrix);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% This section is for orienting the cell so that the antero-posterior pole
% is parallel to the z-axis, with the anterior pole lower than the
% posterior pole on the z-axis.
[x, y, z] = newCoor(x, y, z, antPole, postPole);
cort_x = x(1:num_cort);
cort_y = y(1:num_cort);
cort_z = z(1:num_cort);
oa_x = x(num_cort+1:end);
oa_y = y(num_cort+1:end);
oa_z = z(num_cort+1:end);

antPole = [0 0 norm(antPole - postPole)];
postPole = [0 0 0];

x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if strcmp(slice, 'XYZ')
    K = convhull(x, y, z);
    % plot mesh of convex hull of BB positions
    trimesh(K, x, y, z, 'FaceColor', 'red', 'FaceAlpha', 0.1, ...
        'EdgeColor', 'red', 'EdgeAlpha', 0);
    hold on
    % plot cortical BBs as red filled circles
    scatter3(cort_x, cort_y, cort_z, 80, 'r.');
    % plot oral apparatus BBs as green filled circles
    scatter3(oa_x, oa_y, oa_z, 80, 'g.');
    % plot black straight line between the anterior and posterior poles
    plot3([postPole(1) antPole(1)], [postPole(2) antPole(2)], ...
        [postPole(3) antPole(3)], '-', 'Color', 'black', 'LineWidth', 3);
    
    % color consecutive ciliary rows with different colors for contrast
    colors = {'yellow', 'magenta', 'cyan', 'blue'};
    xy_centroids = zeros(num_rows, 1);
    for i = 1:num_rows
        row = OriginalBBsMatrix(i, :);
        row = row(row ~= 0);
        xy_centroids(i) = atan2(mean(cort_y(row)), mean(cort_x(row)));
    end
    [~, color_order] = sort(xy_centroids, 'ascend');
    for i = 1:num_rows
        row = OriginalBBsMatrix(color_order(i), :);
        row = row(row ~= 0);
        color = colors{mod(i, 4) + 1};
        plot3(cort_x(row), cort_y(row), cort_z(row), 'Color', color, ...
            'LineStyle', '-', 'LineWidth', 1);
    end
    daspect([1 1 (0.125/0.3)]);
    zlabel(slice(3));
elseif strcmp(slice, 'XY')
    K = convhull(x, y);
    fill(x(K), y(K), [1 1 1]);
    hold on
    scatter(cort_x, cort_y, 80, 'r.'); 
    for i = 1:num_rows
        row = OriginalBBsMatrix(i, :);
        row = row(row ~= 0);
        plot(cort_x(row), cort_y(row), 'b-');
    end
elseif strcmp(slice, 'XZ')
    K = convhull(x, z);
    fill(x(K), z(K), [1 1 1]);
    hold on
    scatter(cort_x, cort_z, 80, 'r.');
    for i = 1:num_rows
        row = OriginalBBsMatrix(i, :);
        row = row(row ~= 0);
        plot(cort_x(row), cort_z(row), 'b-');
    end
    daspect([1 (0.125/0.3) 1]);
elseif strcmp(slice, 'YZ')
    K = convhull(y, z);
    fill(y(K), z(K), [1 1 1]);
    hold on
    scatter(cort_y, cort_z, 80, 'r.');
    for i = 1:num_rows
        row = OriginalBBsMatrix(i, :);
        row = row(row ~= 0);
        plot(cort_y(row), cort_z(row), 'b-');
    end
    daspect([1 (0.125/0.3) 1]);
end
title(slice);
xlabel(slice(1));
ylabel(slice(2));
hold off
end

