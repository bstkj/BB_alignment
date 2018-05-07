% -------------------------------------------------------------------------
% [Ben] 03/17/18
% Identifies coordinates for cortical BBs.
% -------------------------------------------------------------------------

function [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, minRow, minCol, maxRow, maxCol] = getBBIdx(I, th)
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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------- OLD VERSION ----------------------------------
% % MIP = max(I, [], 3);
% % MIP_mean = mean(reshape(MIP, [size(MIP,1)*size(MIP,2), 1]));
% % backgroundSubtractedMIP = MIP-MIP_mean;
% % backgroundSubtractedI = I - MIP_mean;
% backgroundSubtractedI = zeros(size(I));
% for i = 1:size(I, 3)
%     background = mean(mean(I(:, :, i)));
%     backgroundSubtractedI(:, :, i) = I(:, :, i) - background;
% end
% 
% % To distinguish between cortical and oral apparatus BBs, a 1um^2 box is
% % centered over the BB maxima and the skew and kurtosis of the intensity
% % histogram are calculatd from a background subtracted image. Cortical BBs
% % have a positive skew and kurtosis, while oral apparatus BBs have a
% % negative skew and kurtosis due to their homogenous background. BBs are
% % excluded if the sum of their skew and kurtosis is less than 1.
% 
% % delete cortical BB
% patchSize = 8; % this used to be 9
% cort = [];
% oa = [];
% for i = 1:length(BB_x_closeToConvexHull)
%     currImg = reshape(backgroundSubtractedI(:, :, BB_z_closeToConvexHull(i)), [height, width]);
%     result = distinguishBBClasses(currImg, BB_x_closeToConvexHull(i), BB_y_closeToConvexHull(i), patchSize);
%     if result == 1
%         cort = [cort i];
%     elseif result == 0
%         oa = [oa i];
%     end
% end
% 
% cort_x = BB_x_closeToConvexHull(cort);
% cort_y = BB_y_closeToConvexHull(cort);
% cort_z = BB_z_closeToConvexHull(cort);
% oa_x = BB_x_closeToConvexHull(oa);
% oa_y = BB_y_closeToConvexHull(oa);
% oa_z = BB_z_closeToConvexHull(oa);
% ---------------------------NEW VERSION ----------------------------------
OA_mask = getPotentialOARegion(I);

[cort, oa] = distinguishBBClasses2(I, BB_x_closeToConvexHull, ...
    BB_y_closeToConvexHull, BB_z_closeToConvexHull, OA_mask);

cort_x = BB_x_closeToConvexHull(cort);
cort_y = BB_y_closeToConvexHull(cort);
cort_z = BB_z_closeToConvexHull(cort);
oa_x = BB_x_closeToConvexHull(oa);
oa_y = BB_y_closeToConvexHull(oa);
oa_z = BB_z_closeToConvexHull(oa);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end