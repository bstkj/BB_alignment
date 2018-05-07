% -------------------------------------------------------------------------
% [Ben] 01/16/18
% Returns the coordinates of the cortical BBs as well as the coordinates of
% the potential OA BBs.  
% -------------------------------------------------------------------------

function [cort_x, cort_y, cort_z, oar_x, oar_y, oar_z] = find_cortBB_OAregion(I, th)

[height, width, ~] = size(I);
[cell, cellArea]=identifyCell(I);
if cellArea < 10000 % equivalent to 10000*(0.125^2) = 156.25 micrometer sq
    [cell, ~] = identifyCell_triangle(I);
end

[row, col] = find(cell == 1);

% boundaries enlarged by 0.125*40 = 5 micrometres
minRow = max(1, min(row) - 40);
minCol = max(1, min(col) - 40);
maxRow = min(height, max(row) + 40);
maxCol = min(width, max(col) + 40);

% cropped image stack
I = I(minRow:maxRow, minCol:maxCol, :);
[height, width, ~] = size(I);

% smoothing by Gaussian kernel
sigma_filter = 1;
I_filtered = planeGaussianFilter(I, sigma_filter);


% get local maxima of 3d image
BW = localMaxima_global(I_filtered); % note that BW is a binary 3D-array

% get potential BBs
[y, x, z] = getPotentialBB_plane(I, BW, th);

k = convhull(x, y, z);
distanceThreshold = 2.25; % unit is micron
minDists = dist2ConvexHull(x, y, z, k);

% 'removed maxima further than 2.25 um from surface of convex hull'
qualified = find(minDists < distanceThreshold);
num_qualified = length(qualified);
BB_x_closeToConvexHull = x(qualified);
BB_y_closeToConvexHull = y(qualified);
BB_z_closeToConvexHull = z(qualified);


MIP = max(I, [], 3);
MIP_mean = mean(reshape(MIP, [size(MIP,1)*size(MIP,2), 1]));
% backgroundSubtractedMIP = MIP-MIP_mean;
backgroundSubtractedI = I-MIP_mean;

% delete cortical BB
patchSize = 9;
potBBsIdx = [];
potOAIdx = [];
for i = 1:num_qualified
    currImg = reshape(backgroundSubtractedI(:, :, BB_z_closeToConvexHull(i)), [height, width]);
    result = distinguishBBClasses(currImg, BB_x_closeToConvexHull(i), BB_y_closeToConvexHull(i), patchSize);
    if result == 1
        potBBsIdx = [potBBsIdx, i];
    else
        potOAIdx = [potOAIdx, i];
    end
end

cort_x = x(potBBsIdx);
cort_y = y(potBBsIdx);
cort_z = z(potBBsIdx);

oar_x = x(potOAIdx);
oar_y = y(potOAIdx);
oar_z = z(potOAIdx);
end

