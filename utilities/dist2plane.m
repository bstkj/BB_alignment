% -------------------------------------------------------------------------
% [Ben] 12/11/17
% Implements the initial step of the BB reassignment described in Jingyi's
% paper on pg 18. Assigns a BB to a ciliary row whose average plane
% (defined as the plane formed by the anterior and posterior poles and the
% centroid of the ciliary row) it is closest to. Returns the distance of
% the BB from that average plane, as well as the corresponding ciliary row
% number.
% -------------------------------------------------------------------------

function [minD, rowIdx] = dist2plane(mtx, x, y, z, ptIdx, antPole, postPole)
numRows = size(mtx, 1);
% numRows is no. of legal ciliary rows detected before BB reassignment
pt2planeDist = zeros(numRows, 1);
pt = [x(ptIdx), y(ptIdx), z(ptIdx)];
for row = 1:numRows
    currRow = mtx(row, :);
    currRow(currRow == 0) = [];
    centroid = [mean(x(currRow)), mean(y(currRow)), mean(z(currRow))];
    % projects pt onto plane defined by anterior and posterior poles, and
    % the centroid of the ciliary row
    pt_proj = projectPtOntoPlane(pt, centroid, antPole, postPole);
    startPts = [centroid; antPole];
    endPts = [pt_proj; postPole];
    pt_intersect = lineIntersect3D(startPts, endPts);
    % What is this condition saying? Maybe to prevent the output from
    % pt2planeDist from blowing up?
    if pt_intersect(1) > min(antPole(1), postPole(1)) && pt_intersect(1) < max(antPole(1), postPole(1))
        pt2planeDist(row) = 1000000;
    else
        pt2planeDist(row) = distance_pt2plane(pt, centroid, antPole, postPole);
    end
end
[minD, rowIdx] = min(pt2planeDist);
end