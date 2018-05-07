function skew_kurt_plot(I, th)
[height, width, ~] = size(I);

[cell, cellArea] = identifyCell(I);
if cellArea < 10000 % equivalent to 10000*(0.125^2) = 156.25 micrometer sq
    [cell, ~] = identifyCell_triangle(I);
end

[row, col] = find(cell == 1);

% boundaries enlarged by 0.125*40 = 5um
% image stack is cropped with the cell mask generated above (enlarged by
% 5um to account for irregularities in the mask)
minRow = max(1, min(row) - 40);
minCol = max(1, min(col) - 40);
maxRow = min(height, max(row) + 40);
maxCol = min(width, max(col) + 40);
I = I(minRow:maxRow, minCol:maxCol, :);
[height, width, ~] = size(I);

% convolved with a low radius (0.12um radius) Gaussian kernel
sigma_filter = 1;
I_filtered = planeGaussianFilter(I, sigma_filter);
% imshow(I_filtered(:, :, 5), [min(min(I_filtered(:, :, 5))), max(max(I_filtered(:, :, 5)))]);
% pause

% I_filtered = imgaussfilt3(I, sigma_filter, 'FilterSize', 0.12/125);

% local maxima are identified using a 0.25um x-y search radius and a 0.6um
% z search radius, which correspond to the approximate dimensions of an
% individual BB
BW = localMaxima_global(I_filtered); % note that BW is a binary 3D-array
% imshow(BW(:, :, 5));
% pause
% BW = localMaxima_global2(I_filtered);
% BW = localMaxima_plane(I_filtered);

% adaptive thresholding approach is used to separate noise from peaks
% corresponding to actual BB maxima
% To calculate the adaptive threshold, the average intensity and standard
% deviation for all peaks is calculated within a rolling average with a z
% range of 1.5um
[y, x, z] = getPotentialBB_plane(I, BW, th);

% idx = find(z == 5);
% imshow(I(:, :, 5), [min(min(I(:, :, 5))), max(max(I(:, :, 5)))]);
% hold on
% scatter(x(idx), y(idx), 80, 'red');
% hold off
% pause

% To remove BB maxima within the interior of the cell, a 3D convex hull is
% generated from the BB maxima and all maxima that are greater than 2.25um
% from the surface of the convex hull are identified and deleted
k = convhull(x, y, z);
distanceThreshold = 2.25; % unit is micron
minDists = dist2ConvexHull(x, y, z, k);
qualified = find(minDists < distanceThreshold);
BB_x_closeToConvexHull = x(qualified);
BB_y_closeToConvexHull = y(qualified);
BB_z_closeToConvexHull = z(qualified);

MIP = max(I, [], 3);
MIP_mean = mean(reshape(MIP, [size(MIP,1)*size(MIP,2), 1]));
% backgroundSubtractedMIP = MIP-MIP_mean;
backgroundSubtractedI = I - MIP_mean;

% To distinguish between cortical and oral apparatus BBs, a 1um^2 box is
% centered over the BB maxima and the skew and kurtosis of the intensity
% histogram are calculatd from a background subtracted image. Cortical BBs
% have a positive skew and kurtosis, while oral apparatus BBs have a
% negative skew and kurtosis due to their homogenous background. BBs are
% excluded if the sum of their skew and kurtosis is less than 1.

% delete cortical BB
patchSize = 10; % this used to be 9
num = length(BB_x_closeToConvexHull);
s = zeros(num, 1);
k = zeros(num, 1);
d = zeros(num, 1);
for i = 1:num
    x_coor = BB_x_closeToConvexHull(i);
    y_coor = BB_y_closeToConvexHull(i);
    z_coor = BB_z_closeToConvexHull(i);
    currImg = reshape(backgroundSubtractedI(:, :, z_coor), [height, width]);
    intensities = currImg((x_coor-patchSize/2):(x_coor+patchSize/2), (y_coor-patchSize/2):(y_coor+patchSize/2));
    intensities_vector = reshape(intensities, [size(intensities, 1)*size(intensities, 2), 1]);
    skew = skewness(intensities_vector);
    kurt = kurtosis(intensities_vector);
    s(i) = skew;
    k(i) = kurt;
    d(i) = z_coor;
end
scatter3(s, k, d, 'r.');
xlabel('skewness');
ylabel('kurtosis');
zlabel('depth');
end