% -------------------------------------------------------------------------
% [Ben] 2/3/18
% Should we merge this with distance_pt2pane.m
% -------------------------------------------------------------------------


function d = dist_pt_to_plane(plane_pt1, plane_pt2, plane_pt3, pt)
% points in input are represented by x,y,z-coor

v1 = plane_pt1 - plane_pt2;
v2 = plane_pt1 - plane_pt3;
v = cross(v1, v2);
d = abs(dot(pt - plane_pt1, v))/norm(v);
end
