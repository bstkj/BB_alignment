% -------------------------------------------------------------------------
% [Ben] 3/16/18 (adapted from Jingyi's getOA.m)
% Gets the coordinates describing the potential OA region in pg 17 of
% Jingyi's paper. These coordinates are stored in the form of a binary
% matrix, and will be subsequently be used in the differentation of
% cortical and oral apparatus BBs.
% -------------------------------------------------------------------------


function OA_mask = getPotentialOARegion(cropped_I)
% 'summed all intensities at each position for all images within image
% stack' pg 17.
SIP = sum(cropped_I, 3); % summed intensity projection
minInt = min(min(SIP));
maxInt = max(max(SIP));
ratio = 0.5;

% 'set a threshold for intensity and a threshold for area' pg 17.
th_int = minInt + ratio*(maxInt - minInt);
th_area = 100;

BW = (SIP > th_int)*1;
[L, n] = bwlabel(BW); % n is # of connected components
nums = zeros(n, 1);

for l=1:n
    % nums stores the 'sizes' of the various connected components
    nums(l) = length(find(L == l));
end

% r stores the indices (relative to nums) of the connected components that
% are 'large' enough
r = find(nums > th_area);
OA_mask = zeros(size(SIP));

for i=1:length(r)
    [row, col] = find(L == r(i));
    % don't correct the indices for cropping here
    for j=1:length(row)
        OA_mask(row(j), col(j)) = 1;
    end
end
end

