imagepath = 'nd2_files/003.nd2';
slice = 'XYZ';
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
% % This section is for orienting the cell so that the antero-posterior pole
% % is parallel to the z-axis, with the anterior pole lower than the
% % posterior pole on the z-axis.
% a = length(cort_x);
% [x, y, z] = newCoor(x, y, z, antPole, postPole);
% cort_x = x(1:a);
% cort_y = y(1:a);
% cort_z = z(1:a);
% oa_x = x(a+1:end);
% oa_y = y(a+1:end);
% oa_z = z(a+1:end);
% [px, py, pz] = newCoor([pole1(1); pole2(1)], [pole1(2); pole2(2)], [pole1(3); pole2(3)], antPole, postPole);
% 
% if pz(1) < pz(2)
%     antPole = [px(2) py(2) pz(2)] - [px(1) py(1) pz(1)];
%     postPole = [0 0 0];
%     cort_x = cort_x - px(1); oa_x = oa_x - px(1);
%     cort_y = cort_y - py(1); oa_y = oa_y - py(1);
%     cort_z = cort_z - pz(1); oa_z = oa_z - pz(1);
% else
%     antPole = [px(1) py(1) pz(1)] - [px(2) py(2) pz(2)];
%     postPole = [0 0 0];
%     cort_x = cort_x - px(2); oa_x = oa_x - px(2);
%     cort_y = cort_y - py(2); oa_y = oa_y - py(2);
%     cort_z = cort_z - pz(2); oa_z = oa_z - pz(2);
% end
% 
% x = vertcat(cort_x, oa_x);
% y = vertcat(cort_y, oa_y);
% z = vertcat(cort_z, oa_z);
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
% % ----------
% % A = horzcat(x, y, z);
% % D = pdist(A);
% % d = max(D);
% % [r, c] = find(squareform(D) == d);
% % antPole = [x(r(1)) y(r(1)) z(r(1))];
% % postPole = [x(c(1)) y(c(1)) z(c(1))];
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
X = horzcat(x, y, z);
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( X, '' );
axes(:, 1) = evecs(:, 1)*radii(1)+center;
axes(:, 4) = -evecs(:, 1)*radii(1)+center;
axes(:, 2) = evecs(:, 2)*radii(2)+center;
axes(:, 5) = -evecs(:, 2)*radii(2)+center;
axes(:, 3) = evecs(:, 3)*radii(3)+center;
axes(:, 6) = -evecs(:, 3)*radii(3)+center;
axes_distance_from_oa = zeros(6, 1);
axes_distance_from_oa(1) = distance_pts(axes(:, 1), OA);
axes_distance_from_oa(4) = distance_pts(axes(:, 4), OA);
axes_distance_from_oa(2) = distance_pts(axes(:, 2), OA);
axes_distance_from_oa(5) = distance_pts(axes(:, 5), OA);
axes_distance_from_oa(3) = distance_pts(axes(:, 3), OA);
axes_distance_from_oa(6) = distance_pts(axes(:, 6), OA);
[~, I] = min(axes_distance_from_oa);
if I <= 3
    antPole = axes(:, I);
    postPole = axes(:, I+3);
else
    antPole = axes(:, I);
    postPole = axes(:, I-3);
end
[x, y, z] = newCoor(x, y, z, antPole, postPole);
a = length(cort_x);
cort_x = x(1:a);
cort_y = y(1:a);
cort_z = z(1:a);
oa_x = x(a+1:end);
oa_y = y(a+1:end);
oa_z = z(a+1:end);
[px, py, pz] = newCoor([antPole(1); postPole(1)], [antPole(2); postPole(2)], ...
    [antPole(3); postPole(3)], antPole, postPole);
% [cx, cy, cz] = newCoor(center(1), center(2), center(3), antPole, postPole);
if pz(1) < pz(2)
    antPole = [px(2) py(2) pz(2)] - [px(1) py(1) pz(1)];
    postPole = [0 0 0];
    cort_x = cort_x - px(1); oa_x = oa_x - px(1);
    cort_y = cort_y - py(1); oa_y = oa_y - py(1);
    cort_z = cort_z - pz(1); oa_z = oa_z - pz(1);
else
    antPole = [px(1) py(1) pz(1)] - [px(2) py(2) pz(2)];
    postPole = [0 0 0];
    cort_x = cort_x - px(2); oa_x = oa_x - px(2);
    cort_y = cort_y - py(2); oa_y = oa_y - py(2);
    cort_z = cort_z - pz(2); oa_z = oa_z - pz(2);
end
x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);
OA = [mean(oa_x), mean(oa_y), mean(oa_z)];

lower_limit = antPole(3)/3;
upper_limit = 2*antPole(3)/3;
% cort_idx = find((cort_z > lower_limit) & (cort_z < upper_limit));
cort_idx = find(cort_z < lower_limit);
cort_x = cort_x(cort_idx);
cort_y = cort_y(cort_idx);
cort_z = cort_z(cort_idx);

if strcmp(slice, 'XYZ')
%     K = convhull(x, y, z);
    % plot mesh of convex hull of BB positions
%     trimesh(K, x, y, z, 'EdgeColor', 'black', 'FaceAlpha', 0.6, ...
%         'FaceColor', 'black');
%     hold on
    % plot cortical BBs as red filled circles 
    scatter3(cort_x, cort_y, cort_z, 160, 'r.');
    hold on
    % plot oral apparatus BBs as green filled circles
    scatter3(oa_x, oa_y, oa_z, 160, 'g.');
    % plot white dotted straight line between the anterior and posterior
    % poles
    plot3([postPole(1), antPole(1)], [postPole(2), antPole(2)], ...
        [postPole(3), antPole(3)], '-', 'Color', [0, 0, 0], 'LineWidth', 2);

    % plot OA centroid
    % scatter3(OA(1), OA(2), OA(3), 2000, 'y.');
    % % plot perpendicular from OA centroid to anterior-posterior axis
    % plot3([OA(1) antPole(1)], [OA(2) antPole(2)], [OA(3) OA(3)], ...
    %     '-', 'Color', [1, 1, 1], 'LineWidth', 2);

    % plots the center of the fitted ellipsoid
    % scatter3(center(1), center(2), center(3), 2000, 'y.');
    % scatter3(cx, cy, cz, 2000, 'y.');
    % 
    % % plots first axis of the fitted ellipsoid
    % axis1 = center+evecs(:, 1)*radii(1);
    % plot3([center(1) axis1(1)], [center(2) axis1(2)], [center(3) axis1(3)], ...
    %     '-', 'Color', [1 0 0], 'LineWidth', 2);
    % 
    % % plots second axis of the fitted ellipsoid
    % axis2 = center+evecs(:, 2)*radii(2);
    % plot3([center(1) axis2(1)], [center(2) axis2(2)], [center(3) axis2(3)], ...
    %     '-', 'Color', [0 1 0], 'LineWidth', 2);
    % 
    % % plots third axis of the fitted ellipsoid
    % axis3 = center+evecs(:, 3)*radii(3);
    % plot3([center(1) axis3(1)], [center(2) axis3(2)], [center(3) axis3(3)], ...
    %     '-', 'Color', [0 0 1], 'LineWidth', 2);

    daspect([1 1 2]); % want data units along x, y, z axes to be equal
    hold off
    zlabel(slice(3));
elseif strcmp(slice, 'XY')
    K = convhull(x, y);
    fill(x(K), y(K), [0.2 0.2 0.2]);
    hold on
    scatter(cort_x, cort_y, 80, 'r.');
    scatter(oa_x, oa_y, 80, 'g.');
%     plot([pole1(1), pole2(1)], [pole1(2), pole2(2)], '--', 'Color', [1, 1, 1]);
    hold off
end