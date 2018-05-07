% -------------------------------------------------------------------------
% [Ben] 12/07/17
% Identifies coordinates for parts of OA region, and returns a vector
% containing the average x, y, z coordinates of the parts identified.
% -------------------------------------------------------------------------

function OA = findOAIdx(I, th)

[height, width, ~] = size(I);
[cell, cellArea] = identifyCell(I);
if cellArea < 10000
    [cell, ~] = identifyCell_triangle(I);
end

[row, col] = find(cell == 1);
minRow = min(row);
minCol = min(col);
maxRow = max(row);
maxCol = max(col);

minRow = max(1, minRow - 40);
minCol = max(1, minCol - 40);
maxRow = min(height, maxRow + 40);
maxCol = min(width, maxCol + 40);


I = I(minRow:maxRow, minCol:maxCol, :);
[height, width, ~] = size(I);

sigma_filter = 1;
I_filtered = planeGaussianFilter(I,sigma_filter);


% get local maxima of 3d image
BW = localMaxima_global(I_filtered);

% get potential BBs
[y, x, z] = getPotentialBB_plane(I, BW, th);

k = convhull(x,y,z);
distanceThreshold = 2.25; % unit is micron
minDists = dist2ConvexHull(x,y,z,k);

qualified = find(minDists < distanceThreshold);
BB_x_closeToConvexHull = x(qualified);
BB_y_closeToConvexHull = y(qualified);
BB_z_closeToConvexHull = z(qualified);


MIP = max(I, [], 3);
MIP_mean = mean(reshape(MIP, [size(MIP,1)*size(MIP,2), 1]));
% backgroundSubtractedMIP = MIP - MIP_mean;
backgroundSubtractedI = I - MIP_mean;

% delete cortical BB
patchSize = 9;
potBBsIdx = [];
for i = 1:length(BB_x_closeToConvexHull)
    currImg = reshape(backgroundSubtractedI(:, :, BB_z_closeToConvexHull(i)), [height,width]);
    result = distinguishBBClasses(currImg, BB_x_closeToConvexHull(i), BB_y_closeToConvexHull(i), patchSize);
    if result == 1
        potBBsIdx = [potBBsIdx, i];
    end
end

allIdx = 1:length(BB_z_closeToConvexHull);
OAIdx = setdiff(allIdx, potBBsIdx); % complement of potBBsIdx

% OA region is defined by the coordinates x(OAIdx), y(OAIdx), z(OAIdx).
% Note that this is not equal to the coordinates of the oral apparatus BBs.
% Rather oral apparatus BBs are identified from this OA region.
OA(1) = round(mean(x(OAIdx)));
OA(2) = round(mean(y(OAIdx)));
OA(3) = round(mean(z(OAIdx)));
end