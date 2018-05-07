% -------------------------------------------------------------------------
% [Ben] 05/01/18 (Written by Ben)
% 3D/2D visualization of basal bodies of a cell for given .nd2 file in
% /nd2_files. The cortical BBs are colored 'red' and the oral apparatus BBs
% are colored 'green'. The line anterior-posterior axis is drawn as a
% dotted white line. This function takes in a path to a .nd2 file, and 
% takes about 13s to run. Depending on what the parameter 'slice' is set
% to, either a 3D mesh plot (of the convex hull of BB positions) or some 2D
% slice of that 3D plot is produced.
% Have not implemented option to control aspect ratio for data.
% -------------------------------------------------------------------------


function visualize_BB(imagepath, slice)
channel = 2; % legacy settings
th = 6; % legacy settings
I = readBioImg(imagepath, channel);

[cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, ~, ~, ~, ~] = getBBIdx(I, th);
x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);

OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);
[pole1, pole2] = findPoles(x, y, z);
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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% This section is for orienting the cell so that the antero-posterior pole
% is parallel to the z-axis, with the anterior pole lower than the
% posterior pole on the z-axis.
num_cort = length(cort_x);
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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% % This section is for rotating the the BB coordinates about the z-axis by
% % some given theta.
% OA = [mean(oa_x), mean(oa_y), mean(oa_z)];
% r = vrrotvec([(OA(1) - antPole(1)), (OA(2) - antPole(2)), 0], [1, 0, 0]);
% theta = r(3)*r(4);
% rotation_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% rotated = rotation_matrix * vertcat(transpose(x), transpose(y));
% x = transpose(rotated(1, :)); y = transpose(rotated(2, :));
% cort_x = x(1:a);
% cort_y = y(1:a);
% oa_x = x(a+1:end);
% oa_y = y(a+1:end);
% OA(1) = mean(oa_x); OA(2) = mean(oa_y);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if strcmp(slice, 'XYZ')
    K = convhull(x, y, z);
    % plot mesh of convex hull of BB positions
    trimesh(K, x, y, z, 'FaceColor', 'blue', 'FaceAlpha', 0.2, ...
        'EdgeColor', 'blue', 'EdgeAlpha', 0.4);
    hold on
    % plot cortical BBs as red filled circles 
    scatter3(cort_x, cort_y, cort_z, 160, 'r.');
    % plot oral apparatus BBs as green filled circles
    scatter3(oa_x, oa_y, oa_z, 160, 'g.');
    % plot yellow straight line between the anterior and posterior poles
    plot3([postPole(1), antPole(1)], [postPole(2), antPole(2)], ...
        [postPole(3), antPole(3)], '-', 'Color', 'yellow', 'LineWidth', 3);
%     % plot OA centroid
%     scatter3(OA(1), OA(2), OA(3), 2000, 'y.');
%     % plot perpendicular from OA centroid to anterior-posterior axis
%     plot3([OA(1) antPole(1)], [OA(2) antPole(2)], [OA(3) OA(3)], ...
%         '-', 'Color', [1, 1, 1], 'LineWidth', 2);
    daspect([1 1 (0.125/0.3)]);
    hold off
    zlabel(slice(3));
elseif strcmp(slice, 'XY')
    K = convhull(x, y);
    fill(x(K), y(K), [0.5 0.5 1]);
    hold on
    scatter(cort_x, cort_y, 80, 'r.');
    scatter(oa_x, oa_y, 80, 'g.');
    hold off
elseif strcmp(slice, 'XZ')
    K = convhull(x, z);
    fill(x(K), z(K), [0.5 0.5 1]);
    hold on
    scatter(cort_x, cort_z, 80, 'r.');
    scatter(oa_x, oa_z, 80, 'g.');
    daspect([1 (0.125/0.3) 1]);
    hold off
elseif strcmp(slice, 'YZ')
    K = convhull(y, z);
    fill(y(K), z(K), [0.5 0.5 1]);
    hold on
    scatter(cort_y, cort_z, 80, 'r.');
    scatter(oa_y, oa_z, 80, 'g.');
    daspect([1 (0.125/0.3) 1]);
    hold off
end
title(slice);
xlabel(slice(1));
ylabel(slice(2));
end

