% -------------------------------------------------------------------------
% [Ben] 1/31/18
% Calculates perpendicular distance from pt to plane defined by plane_pt1,
% plane_pt2, and plane_pt3. 
% -------------------------------------------------------------------------

function d = distance_pt2plane(pt, plane_pt1, plane_pt2, plane_pt3)
% why is the scaling different? e.g: 0.1278 instead of 0.125
% scale_vector = [0.1278, 0.1278, 0.3];
scale_vector = [0.125 0.125 0.3];

v1 = (plane_pt1 - plane_pt2) .* scale_vector;
v2 = (plane_pt1 - plane_pt3) .* scale_vector;
% normal vector to the plane
n = cross(v1, v2)/norm(cross(v1, v2)); 
% absolute value of dot product is the perpendicular distance to the plane
d = abs(sum(((pt - plane_pt1) .* scale_vector) .* n));
end



