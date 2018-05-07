% -------------------------------------------------------------------------
% [Ben] (adapted from Jingyi's code by Ben) 01/02/18
% Performs the cell detection outlined in pg 16 of Jinyi's paper.
% Extracting MIP, applying several filters, binarizing image, searching for
% largest object in resultant image satisfying several conditions.
% Returns mask of detected cell, together with area of mask.
% -------------------------------------------------------------------------

function [cell, cellArea] = identifyCell(oriI)

% MIP created and mean pixel intensity is subtracted from MIP
[height, width, ~] = size(oriI);
oriMIP = max(oriI, [], 3); % return maxes along dim 3 (depth)
totalMax = max(max(oriMIP));
totalMin = min(min(oriMIP));
% For grayscale images: For double arrays, values range from [0, 1]. For
% uint8, values range from [0, 255].

% MIP = uint8((oriMIP-totalMin)/(totalMax-totalMin)*255);
% MIP_mean = sum(sum(MIP))/(size(MIP, 1)*size(MIP, 2));
MIP = uint8(zeros(height, width));
for i = 1:height
    for j = 1:width
        MIP(i, j) = (oriMIP(i, j)-totalMin)/(totalMax - totalMin)*255;
    end
end
MIP_mean = mean(reshape(MIP, [size(MIP, 1)*size(MIP, 2), 1]));
backgroundSubtractedMIP = MIP - MIP_mean;

% Background subtracted MIP is convolved with large radius Gaussian kernel
% (1um radius) to get larger cellular scale features
B = imgaussfilt(backgroundSubtractedMIP, 1); % standard deviation of 1
% B = imgaussfilt(backgroundSubtractedMIP, sigma, 'FilterSize', 1/0.125);


% Large features are separated by convolving with small radius Laplacian
% (0.12um radius) 
% and further smoothed with another large radius Gaussian kernel
h = fspecial('log'); % Laplacian of Gaussian filter
% h = fspecial('log', 0.12/0.125, sigma); 
C = imfilter(B, h);
D = imgaussfilt(C, 8);
% D = imgaussfilt(C, sigma, 'FilterSize', 1/0.125);

% resulting MIP is thresholded using Otsu's method and segmented objects
% are filtered based on their shape
level = graythresh(D); 
BW = imbinarize(D, level);

L = bwlabel(BW); % label connected components
areaSize = zeros(max(max(L)), 1);
for l = 1:max(max(L))
    [y, x] = find(L == l);
    try
        K = convhull(x,y);
        areaSize(l)=polyarea(x(K), y(K));
        
        % area represented by 1 pixel seems to be (0.125^2) 
        % objects > 2000um^2 in 2D area
        if areaSize(l) > 2000/(0.125^2)
            areaSize(l) = 0;
        end

        % feret diameter > 60um
        currImg = zeros(height, width);
        for i = 1:length(x)
            currImg(y(i), x(i)) = 1;
        end
        
        currImg = imfill(currImg);
        [fd, ~] = imFeretDiameter(currImg);
        if fd > 60/0.125
            areaSize(l) = 0;
        end
        
        % circularity < 0.85
        stats = regionprops(currImg, 'Perimeter', 'Area');   
  
        circularity = 4*pi*stats.Area/stats.Perimeter^2;
        if circularity < 0.85 % perfect circle has circularity == 1
            areaSize(l) = 0;
        end 
    catch
        areaSize(l) = 0;
    end
end

% largest remaining object that is retained is the cell of interest
cellArea = max(areaSize);
label = find(areaSize == cellArea, 1);
[y, x] = find(L == label);
cell = zeros(height, width);

for i=1:length(x)
    cell(y(i), x(i)) = 1;
end

cell = imfill(cell);
end