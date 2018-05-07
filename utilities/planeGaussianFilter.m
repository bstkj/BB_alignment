% -------------------------------------------------------------------------
% [Ben] 2/2/18
% Applies Gaussian filter to each xy plane of the 3D image (pixel) array.
% Each 3D array is formed by 2D (xy plane) matrices stacked on top of each
% other.
% Should we be applying a 2D Gaussian filter to each xy plane of the 3D
% image, or should we be applying a 3D Gaussian filter instead??
% -------------------------------------------------------------------------
function filter=planeGaussianFilter(I, sigma)
[m, n, l]=size(I);
filter = zeros(m, n, l);
for i = 1:l
    currPlane = reshape(I(:, :, i), [m, n]);
%     currFilter = imgaussfilt(currPlane, sigma);
    currFilter = imgaussfilt(currPlane, sigma, 'FilterSize', 3);
    filter(:, :, i) = currFilter;
end
end