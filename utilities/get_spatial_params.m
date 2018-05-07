% -------------------------------------------------------------------------
% [Ben] 05/03/18
% For the time being this only calculates the number of BBs in each of the
% 16 regions that we divide the cell into.
% -------------------------------------------------------------------------

function region_count = get_spatial_params(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, antPole)

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Rotates the BB coordinates about the z-axis such that the perpendicular
% from the OA centroid to the antero-posterior axis is aligned with the
% positive x-axis
OA = [mean(oa_x), mean(oa_y), mean(oa_z)];
r = vrrotvec([(OA(1) - antPole(1)), (OA(2) - antPole(2)), 0], [1, 0, 0]);
theta = r(3)*r(4);
rotation_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
rotated = rotation_matrix * vertcat(transpose(cort_x), transpose(cort_y));

cort_x = transpose(rotated(1, :)); 
cort_y = transpose(rotated(2, :));
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

plane1_z = 0.25*antPole(3);
plane2_z = 0.5*antPole(3);
plane3_z = 0.75*antPole(3);

% regions start from the quadrants below plane1, in the order of quadrant
% 1, 2, 3, 4 and then go to the quadrants between plane2 and plane1 ... and
% so on.
region_count = zeros(1, 16);
% section no. (section 1 is under plane1, section 2 is between plane 2 and
% plane 3 ... and so on)
% quad no. (quad 1, 2, 3, 4 are defined according to common trig
% conventions)

section1_idx = find(cort_z < plane1_z);

section2_idx = find(cort_z < plane2_z); 
section2_idx = setdiff(section2_idx, section1_idx);

section3_idx = find(cort_z < plane3_z); 
section3_idx = setdiff(setdiff(section3_idx, section2_idx), section1_idx);

section4_idx = find(cort_z >= plane3_z);

for i = 1:4
    if i == 1
        idx = section1_idx;
    elseif i == 2
        idx = section2_idx;
    elseif i == 3
        idx = section3_idx;
    else
        idx = section4_idx;
    end
    k = 4*(i-1);
    for j = 1:length(idx)
        xj = cort_x(idx(j)); yj = cort_y(idx(j));
        if xj >= 0
            if yj >= 0
                region_count(1+k) = region_count(1+k) + 1;
            else
                region_count(4+k) = region_count(4+k) + 1;
            end
        else
            if yj >= 0
                region_count(2+k) = region_count(2+k) + 1;
            else
                region_count(3+k) = region_count(3+k) + 1;
            end
        end
    end
end

end

