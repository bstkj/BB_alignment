% -------------------------------------------------------------------------
% [Ben] 12/07/17
% Calculates the minimum distance between each potential BB point, and the
% convex hull created from the set of potential BB points. Returns the
% minimum distances calculated. 
% Note that by 'potential BB point', I am referring to the output of
% 'getPotentialBB_plane'. 
% -------------------------------------------------------------------------

function minDists = dist2ConvexHull(x, y, z, k)

x = x*0.125;
y = y*0.125;
z = z*0.3;

minDists = zeros(length(x), 1);
for BB = 1:length(x)
    pt = [x(BB), y(BB), z(BB)];
    currDists = zeros(size(k, 1), 1);
    for i = 1:size(k, 1)
        pt1 = [x(k(i, 1)), y(k(i, 1)), z(k(i, 1))];
        pt2 = [x(k(i, 2)), y(k(i, 2)), z(k(i, 2))];
        pt3 = [x(k(i, 3)), y(k(i, 3)), z(k(i, 3))];
        d = dist_pt_to_plane(pt1, pt2, pt3, pt);
        currDists(i) = d;
    end
    minDists(BB) = min(currDists);
end
end