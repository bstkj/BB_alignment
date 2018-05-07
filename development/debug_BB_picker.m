% -------------------------------------------------------------------------
% [Ben] 4/1/18
% Still thinking of the range of functionalities this should include. To
% start with, should have the ability to display the different stages
% through which BBs are eventually picked out. 
% -------------------------------------------------------------------------

function debug_BB_picker(imagepath, zplane)
channel = 2; % legacy settings
th = 6; % legacy settings
I = readBioImg(imagepath, channel);

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


MIP = max(I, [], 3);

% convolved with a low radius (0.12um radius) Gaussian kernel
sigma_filter = 1;
I_filtered = planeGaussianFilter(I, sigma_filter);

imshow(I_filtered(:, :, zplane), [min(min(I_filtered(:, :, zplane))) ...
    max(max(I_filtered(:, :, zplane)))])
pause

background1 = imopen(I_filtered, strel('sphere', 3));

I2 = I_filtered - background1;
imshow(I2(:, :, zplane), [min(min(I2(:, :, zplane))) max(max(I2(:, :, zplane)))]);
pause

J2 = MIP - background1(:, :, zplane);
imshow(J2, [min(min(J2)) max(max(J2))]);
pause

imshow(background1(:, :, zplane), [min(min(background1(:, :, zplane))) ...
    max(max(background1(:, :, zplane)))]);
pause

BW = localMaxima_global(I_filtered); % note that BW is a binary 3D-array

% adaptive thresholding approach is used to separate noise from peaks
% corresponding to actual BB maxima
% To calculate the adaptive threshold, the average intensity and standard
% deviation for all peaks is calculated within a rolling average with a z
% range of 1.5um
[y, x, z] = getPotentialBB_plane(I2, BW, th);

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
imshow(I(:, :, zplane), [min(min(I(:, :, zplane))) ...
    max(max(I(:, :, zplane)))])
hold on
scatter3(x(idx), y(idx), z(idx), 'red');
hold off
pause
end

