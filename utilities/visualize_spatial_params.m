% -------------------------------------------------------------------------
% [Ben] 4/12/18
% For the time being this creates a heatmap from the output of
% get_spatial_params.m
% -------------------------------------------------------------------------


function visualize_spatial_params(region_count)
region_count_mat = [region_count(1:4); region_count(5:8); ...
    region_count(9:12); region_count(13:16)];
colormap('hot');
image(region_count_mat);
xlabel('Quadrant no.', 'FontSize', 15); ylabel('Section no.', 'FontSize', 15);
xticks([1 2 3 4]); yticks([1 2 3 4]);
title('Regional BB counts', 'FontSize', 15);
c = colorbar;
c.Label.String = 'no. of BBs';
c.Label.FontSize = 15;
end

