% -------------------------------------------------------------------------
% [Ben] 04/24/18
% Produces graphs comparing the 'initially-assigned' ciliary rows based on
% statistics such as no. of BBs per row, ave. distance between consecutive
% BBs for each row, ave. angular deviation between consecutive BBs for each
% row ...
% -------------------------------------------------------------------------

function visualize_initial_CR_stats(imagepath, load_saved, matfile_path)
channel = 2; % legacy settings
th = 6; % legacy settings


I = readBioImg(imagepath, channel);
minBBsInRow = 10;
[cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, minRow, minCol, maxRow, maxCol] = getBBIdx(I, th);
x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);
OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);
cropped_I = I(minRow:maxRow, minCol:maxCol, :);

OA_mask = getPotentialOARegion(cropped_I);
[oa_mask_y, oa_mask_x] = find(OA_mask);
oa_mask_z = repelem(10, size(oa_mask_y, 1), 1);

if load_saved
    load(matfile_path);
    updatedTraceback = initialMatrix;
else
    [~, updatedTraceback] = getInitialCiliaryRows(x, y, z, OA, minBBsInRow);
end

% histogram comparing no. of BBs in each row
histogram(sum(updatedTraceback > 0, 2), 30, 'FaceColor', 'red');
label = sprintf('Total no. of basal bodies: %d', sum(sum(updatedTraceback > 0)));
legend({label}, 'FontSize', 12);
xlabel('No. of basal bodies');
ylabel('No. of ciliary rows');
pause

% histogram comparing the ave. distance between consecutive BBs for each
% row
ave_distances = zeros(size(updatedTraceback, 1), 1);
for i = 1:size(updatedTraceback, 1)
    row = updatedTraceback(i, :);
    row = row(row ~= 0);
    distances = zeros(length(row)-1, 1);
    for j = 1:(length(row)-1)
        bb_1 = [x(row(j)) y(row(j)) z(row(j))];
        bb_2 = [x(row(j+1)) y(row(j+1)) z(row(j+1))];
        distances(j) = distance_pts(bb_1, bb_2);
    end
    ave_distances(i) = mean(distances);
end
histogram(ave_distances, 30, 'FaceColor', 'green');
label = sprintf('Total no. of ciliary rows: %d', size(updatedTraceback, 1));
legend({label}, 'FontSize', 12);
xlabel('Ave. distance between consecutive basal bodies (um)');
ylabel('No. of ciliary rows');
pause

% histogram comparing the ave. angular deviation between consecutive BBs
% for each row
ave_angular_devs = zeros(size(updatedTraceback, 1), 1);
for i = 1:size(updatedTraceback, 1)
    row = updatedTraceback(i, :);
    row = row(row ~= 0);
    displacement_vectors = cell(length(row)-1, 1);
    for j = 1:(length(row)-1)
        bb_1 = [x(row(j)) y(row(j)) z(row(j))];
        bb_2 = [x(row(j+1)) y(row(j+1)) z(row(j+1))];
        displacement_vectors{j} = bb_2 - bb_1;
    end
    angular_devs = zeros(length(displacement_vectors)-1, 1);
    for j = 1:(length(displacement_vectors)-1)
        v_1 = displacement_vectors{j};
        v_2 = displacement_vectors{j+1};
        if norm(cross(v_1, v_2)) % v_1 and v_2 are not parallel
            cos_theta = dot(v_1, v_2)/(norm(v_1)*norm(v_2));
            theta = acosd(cos_theta);
        else % this conditional was added to take care of precision issues
            theta = 0;
        end
        angular_devs(j) = theta;
    end
    ave_angular_devs(i) = mean(angular_devs);
end
histogram(ave_angular_devs, 30, 'FaceColor', 'blue');
xlabel('Ave. angular deviation between consecutive basal bodies (deg)');
ylabel('No. of ciliary rows');
pause

% scatter plot of ave. angular deviation between consecutive basal bodies
% vs. no. of basal bodies in given row
scatter(sum(updatedTraceback > 0, 2), ave_angular_devs, 200, 'r.');
xlabel('No. of basal bodies');
ylabel('Ave. angular deviation between consecutive basal bodies (deg)');
pause

% 3d scatter plot of ave. angular deviation between consecutive basal
% bodies vs. no. of basal bodies in given row vs. ave. distance between
% consecutive basal bodies in given row
scatter3(sum(updatedTraceback > 0, 2), ave_angular_devs, ave_distances, 'filled', 'green');
xlabel('No. of basal bodies');
ylabel('Ave. angular deviation between consecutive basal bodies (deg)');
zlabel('Ave. distance between consecutive basal bodies (um)');
pause

% reorient the cell and produce heat maps + density plots
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

[rot_x, rot_y, rot_z] = newCoor(x, y, z, antPole, postPole);

[~, convexhull_volume] = convhull(rot_x, rot_y, rot_z);
convexhull_volume = convexhull_volume*(0.125*0.125*0.3);
fprintf('Volume bounded by convex hull: %d (um^3)\n', convexhull_volume);

[rot_oa_mask_x, rot_oa_mask_y, rot_oa_mask_z] = newCoor(oa_mask_x, ...
    oa_mask_y, oa_mask_z, antPole, postPole);
scatter3(rot_x, rot_y, rot_z, 80, 'r.');
hold on
scatter3(rot_oa_mask_x, rot_oa_mask_y, rot_oa_mask_z, 'b.');
for i = 1:size(updatedTraceback, 1)
    row = updatedTraceback(i, :);
    row = row(row ~= 0);
    plot3(rot_x(row), rot_y(row), rot_z(row), 'c-');
end

% calculate length/width/volume
X = horzcat(rot_x, rot_y, rot_z);
[center, radii, evecs, ~, ~] = ellipsoid_fit( X, '' );

estimated_volume = (4/3)*pi*radii(1)*radii(2)*radii(3);
estimated_volume = estimated_volume*(0.125*0.125*0.3);
fprintf('Volume bounded by fitted ellipsoid: %d (um^3)\n', estimated_volume);

%     plots the center of the fitted ellipsoid
scatter3(center(1), center(2), center(3), 2000, 'y.');

% plots first axis of the fitted ellipsoid
axis1 = center+evecs(:, 1)*radii(1);
plot3([center(1) axis1(1)], [center(2) axis1(2)], [center(3) axis1(3)], ...
    '-', 'Color', [1 0 0], 'LineWidth', 2);

% plots second axis of the fitted ellipsoid
axis2 = center+evecs(:, 2)*radii(2);
plot3([center(1) axis2(1)], [center(2) axis2(2)], [center(3) axis2(3)], ...
    '-', 'Color', [0 1 0], 'LineWidth', 2);

% plots third axis of the fitted ellipsoid
axis3 = center+evecs(:, 3)*radii(3);
plot3([center(1) axis3(1)], [center(2) axis3(2)], [center(3) axis3(3)], ...
    '-', 'Color', [0 0 1], 'LineWidth', 2);

hold off
pause

rotation1 = vrrotvec(evecs(:, 3), [0 0 1]);
rot_mat1 = vrrotvec2mat(rotation1);
% rotations in 3D are generally not commutative
rotation2 = vrrotvec(rot_mat1*evecs(:, 1), [1 0 0]);
rot_mat2 = vrrotvec2mat(rotation2);
X2 = (rot_mat2*rot_mat1*vertcat(rot_x', rot_y', rot_z'));
X2 = X2 - (rot_mat2*rot_mat1*center);
scatter3(X2(1, :), X2(2, :), X2(3, :), 80, 'r.');
hold on
[e_x, e_y, e_z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 30);
surf(e_x, e_y, e_z, 'FaceColor', 'green', 'FaceAlpha', 0.2, ...
    'EdgeColor', 'green', 'EdgeAlpha', 0.4);

K = convhull(X2(1, :), X2(2, :), X2(3, :));
trimesh(K, X2(1, :), X2(2, :), X2(3, :), 'FaceColor', 'blue', ...
    'FaceAlpha', 0.2, ... 
    'EdgeColor', 'blue', 'EdgeAlpha', 0.4);

hold off
pause

% spatial distribution of BBs

% region_count = get_spatial_params(imagepath);
% visualize_spatial_params(region_count);
% pause

% polarity stuff 

% publish various figures into a report

end

