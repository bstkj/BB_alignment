% -------------------------------------------------------------------------
% [Ben] 03/17/18 (largely adapted from getBBIdx.m)
% Superposes XY slice of basal body positions onto maximum intensity
% projection for given .nd2 file in /nd2_files. The cortical BBs are
% colored 'red' while the oral apparatus BBs are colored 'green'.
% show_cropped is a boolean parameter. If set to true, then the resultant
% image will be cropped to a box just bigger than, and containing, the
% cell. show_OA_region is a boolean parameter. If set to true, then the
% resultant image will also show all the potential OA BBs in a different
% color from the cortical BBs. Takes about 13s to run.
% For example: superpose_BB_MIP('nd2_files/001.nd2', true)
% -------------------------------------------------------------------------


function superpose_BB_MIP(imagepath, show_cropped)
channel = 2; % legacy settings
th = 6; % legacy settings
I = readBioImg(imagepath, channel);

[cort_x, cort_y, ~, oa_x, oa_y, ~, minRow, minCol, maxRow, maxCol] = getBBIdx(I, th);

MIP = max(I, [], 3); 
totalMax = max(max(MIP));
totalMin = min(min(MIP));
norm_MIP = (MIP-totalMin)/(totalMax-totalMin);

if show_cropped
    imshow(norm_MIP(minRow:maxRow, minCol:maxCol));
    hold on
    scatter(cort_x, cort_y, 'r.');
    scatter(oa_x, oa_y, 'g.');
    hold off
else
   imshow(norm_MIP);
   hold on
   scatter(cort_x+minCol, cort_y+minRow, 'r.');
   scatter(oa_x+minCol, oa_y+minRow, 'g.');
   hold off
end
end

