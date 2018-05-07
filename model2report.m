% -------------------------------------------------------------------------
% [Ben]
% Generates figures pertaining to the attributes listed below, saves them
% as .png files 

% no. of BBs
% no. of ciliary rows
% no. of BBs per row
% distribution of ave. distance between consecutive BBs for each row
% distribution of ave. angular deviation between consecutive BBs for each
% row
% heat map of no. of BBs in different regions of the cell
% length/width/volume if we fit an ellipsoid to the 'framework' formed by
% the BBs
% polarity marker stuff
% -------------------------------------------------------------------------

function model2report(modelpath)
% should load model struct with the following attributes: cort_x,
% cort_y, cort_z, oa_x, oa_y, oa_z, antPole, postPole, final_traceback
data = load(modelpath);
% these coordinates should already have been rotated and translated as
% in visualize_initial_CR
cort_x = data.model.cort_x; oa_x = data.model.oa_x;
cort_y = data.model.cort_y; oa_y = data.model.oa_y;
cort_z = data.model.cort_z; oa_z = data.model.oa_z;
antPole = data.model.antPole;
final_traceback = data.model.final_traceback;


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% prints total no. of cortical BBs picked out to screen
fprintf('Total no. of cortical basal bodies: %d\n', ...
    sum(sum(final_traceback > 0)));
% prints total no. of ciliary rows to screen
fprintf('Total no. of ciliary rows: %d\n', size(final_traceback, 1));
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% histogram comparing no. of BBs in each row
histogram(sum(final_traceback > 0, 2), 30, 'FaceColor', 'red');
xlabel('No. of basal bodies');
ylabel('No. of ciliary rows');
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% distribution of ave. distance between consecutive BBs for each row
ave_distances = zeros(size(final_traceback, 1), 1);
for i = 1:size(final_traceback, 1)
    row = final_traceback(i, :);
    row = row(row ~= 0);
    distances = zeros(length(row)-1, 1);
    for j = 1:(length(row)-1)
        bb_1 = [cort_x(row(j)) cort_y(row(j)) cort_z(row(j))];
        bb_2 = [cort_x(row(j+1)) cort_y(row(j+1)) cort_z(row(j+1))];
        distances(j) = distance_pts(bb_1, bb_2);
    end
    ave_distances(i) = mean(distances);
end

histogram(ave_distances, 30, 'FaceColor', 'green');
xlabel('Ave. distance between consecutive basal bodies (um)');
ylabel('No. of ciliary rows');
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% histogram comparing the ave. angular deviation between consecutive BBs
% for each row
ave_angular_devs = zeros(size(final_traceback, 1), 1);
for i = 1:size(final_traceback, 1)
    row = final_traceback(i, :);
    row = row(row ~= 0);
    displacement_vectors = cell(length(row)-1, 1);
    for j = 1:(length(row)-1)
        bb_1 = [cort_x(row(j)) cort_y(row(j)) cort_z(row(j))];
        bb_2 = [cort_x(row(j+1)) cort_y(row(j+1)) cort_z(row(j+1))];
        displacement_vectors{j} = bb_2 - bb_1;
    end
    angular_devs = zeros(length(displacement_vectors)-1, 1);
    for j = 1:(length(displacement_vectors)-1)
        v_1 = displacement_vectors{j};
        v_2 = displacement_vectors{j+1};
        cos_theta = dot(v_1, v_2)/(norm(v_1)*norm(v_2));
        theta = acosd(cos_theta);
%         if norm(cross(v_1, v_2)) % v_1 and v_2 are not parallel
%             cos_theta = dot(v_1, v_2)/(norm(v_1)*norm(v_2));
%             theta = acosd(cos_theta);
%         else % this conditional was added to take care of precision issues
%             theta = 0;
%         end
        angular_devs(j) = theta;
    end
    ave_angular_devs(i) = mean(angular_devs);
end
% the line below was added as a fix to some precision issues that were
% producing imaginary parts to the thetas calculated
ave_angular_devs = real(ave_angular_devs);
histogram(ave_angular_devs , 30, 'FaceColor', 'blue');
xlabel('Ave. angular deviation between consecutive basal bodies (deg)');
ylabel('No. of ciliary rows');
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% heat map of no. of BBs in different regions of the cell (spatial
% distribution of the BBs)
region_count = get_spatial_params(cort_x, cort_y, cort_z, oa_x, oa_y, ...
    oa_z, antPole);
visualize_spatial_params(region_count);
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% length/width/volume if we fit the convex hull / an ellipsoid to the
% 'framework' formed by the BBs

x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);

X = horzcat(x, y, z);
[center, radii, evecs, ~, ~] = ellipsoid_fit( X, '' );

% not sure how the columns of evecs are arranged, or if they are even
% sorted in any way
rotation1 = vrrotvec(evecs(:, 3), [0 0 1]);
rot_mat1 = vrrotvec2mat(rotation1);
% rotations in 3D are generally not commutative
rotation2 = vrrotvec(rot_mat1*evecs(:, 1), [1 0 0]);
rot_mat2 = vrrotvec2mat(rotation2);
X2 = (rot_mat2*rot_mat1*vertcat(x', y', z'));
X2 = X2 - (rot_mat2*rot_mat1*center);
scatter3(X2(1, :), X2(2, :), X2(3, :), 80, 'r.');
hold on
[e_x, e_y, e_z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 30);
surf(e_x, e_y, e_z, 'FaceColor', 'green', 'FaceAlpha', 0.2, ...
    'EdgeColor', 'green', 'EdgeAlpha', 0.4);

[K, convexhull_volume] = convhull(X2(1, :), X2(2, :), X2(3, :));
trimesh(K, X2(1, :), X2(2, :), X2(3, :), 'FaceColor', 'blue', ...
    'FaceAlpha', 0.2, ... 
    'EdgeColor', 'blue', 'EdgeAlpha', 0.4);
legend({'All BBs (cort & oa)', 'Ellipsoid fit', 'Convex hull'}, 'FontSize', 12);
hold off

% Multiplying by (0.125*0.125*0.3) is an estimation, but gives a number in
% the right ball park
convexhull_volume = convexhull_volume*(0.125*0.125*0.3);
fprintf('Volume bounded by convex hull: %d (um^3)\n', convexhull_volume);

% Volume of ellipsoid
estimated_volume = (4/3)*pi*radii(1)*radii(2)*radii(3);
estimated_volume = estimated_volume*(0.125*0.125*0.3);
fprintf('Volume bounded by fitted ellipsoid: %d (um^3)\n', estimated_volume);

% Dimensions of the fitted ellipsoid
fprintf('Lengths of principal axes of the fitted ellipsoid (um): %d, %d, %d', ...
    radii(1), radii(2), radii(3));
pause
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< <<<<<<

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% polarity marker stuff

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end

