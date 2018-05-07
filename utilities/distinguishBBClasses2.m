% -------------------------------------------------------------------------
% [Ben] 3/16/18
% Returns a vector of indices (relative to the vector x/y/z) corresponding
% to cortical BBs, and a similar vector corresponding to oral apparatus
% BBs. Uses the approach described in pg 17-18 of Jingyi's paper, which
% involves having previously identified the general location/boundaries of
% the oral apparatus (this information is passed in as OA_mask). What is
% then done is to use the characteristics of the 'average' cortical BB to
% filter the potential oral apparatus BBs. This is in contrast. 
% -------------------------------------------------------------------------


function [cort, oa] = distinguishBBClasses2(img, x, y, z, OA_mask)
% img and OA_mask have the same dimensions

r = 1.01;
cort = [];
oa = [];
depth = size(img, 3);

for i=1:depth
    idx = find(z == i);
    BB_x = x(idx);
    BB_y = y(idx);
    % 'For each plane, we calculate the average intensity of all pixels
    % within a 0.25um*0.25um*0.6um box centered over BB candidates located
    % at a 0-value position in OA_mask' pg 17.
    cort_intensities = [];
    cort_idx = [];
    pot_oa_intensities = [];
    pot_oa_idx = [];
    for j=1:length(idx)
        % instead of a 2 x 2 x 2 box, we use a 3 x 3 x 3 box
        box = img((BB_y(j)-1):(BB_y(j)+1), (BB_x(j)-1):(BB_x(j)+1), max(1,i-1):min(depth,i+1));
        [x_len, y_len, z_len] = size(box);
        ave_intensity = sum(sum(sum(box)))/(x_len*y_len*z_len);
        
        if OA_mask(BB_y(j), BB_x(j)) == 0
            cort_intensities = [cort_intensities ave_intensity];
            cort_idx = [cort_idx idx(j)];
        else
            pot_oa_intensities = [pot_oa_intensities ave_intensity];
            pot_oa_idx = [pot_oa_idx idx(j)]; % indices relative to the vector BB_x
        end
    end
    % 'For BBs detected within the OA region, if the average intensity of
    % its corresponding box exceeds the average intensity of cortical BBs
    % in the given plane times a ratio, r, we treat it as an oral apparatus
    % BB' pg 18.
    ave_cort_intensity = mean(cort_intensities);
    oa_idx = pot_oa_idx((pot_oa_intensities/ave_cort_intensity) > r);
    
    cort = [cort cort_idx];
    oa = [oa oa_idx];
end
end

