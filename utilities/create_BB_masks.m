% -------------------------------------------------------------------------
% [Ben] 05/02/18
% Creates a mask indicating the locations of the detected cortical and oral
% apparatus BBs at every z-slice of I. This function is called by
% create_overlay.m. Note that cortical BBs are marked with 1's, whereas
% oral apparatus BBs are marked with 2's.
% -------------------------------------------------------------------------

function masks = create_BB_masks(I, th)
[cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, minRow, minCol, ~, ~] = getBBIdx(I, th);
[rows, cols, layers] = size(I);
masks = zeros(rows, cols, layers);
% cortical BBs
for i = 1:length(cort_x)
    masks(cort_y(i)+minRow, cort_x(i)+minCol, cort_z(i)) = 1;
end
% oral apparatus BBs
for i=1:length(oa_x)
    masks(oa_y(i)+minRow, oa_x(i)+minCol, oa_z(i)) = 2;
end
end

