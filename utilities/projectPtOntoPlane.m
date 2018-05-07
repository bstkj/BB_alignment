% -------------------------------------------------------------------------
% [Ben] 2/2/18
% Returns the projection of pt onto the plane defined by pt1, pt2, and pt3.
% pt is projected to the point in the plane that has the least distance
% from it.
% -------------------------------------------------------------------------

function [projectionPt] = projectPtOntoPlane(pt, pt1, pt2, pt3)
v1 = pt1 - pt2;
v2 = pt1 - pt3;
n = cross(v1, v2)/norm(cross(v1, v2));
projectionPt = pt - dot(pt - pt1, n)*n;
end