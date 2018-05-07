% -------------------------------------------------------------------------
% [Ben] 12/19/17 (adapted from Jingyi's code, by Ben)
% Returns the optimized ciliary row configuration, as well as the
% corresponding F-score
% -------------------------------------------------------------------------

function [f_final, initialMatrix] = getCiliaryRows(x, y, z, antPole, postPole, minBBsInRow)
% To identify the poles of the cell, the maximum distance between each BB
% and every other BB is determined and a list of the 10 greatest distances
% is created. These distances correspond to the 10 most anterior BBs and
% the 10 most posterior BBs and the centroid of each of these clusters of
% BBs is taken to be the anterior and posterior pole respectively.

num_BBs = length(x);
dist2Ant = zeros(num_BBs, 1);
dist2Post = zeros(num_BBs, 1);
for i = 1:num_BBs
    BB = [x(i), y(i), z(i)];
    dist2Ant(i) = distance_pts(BB, antPole);
    dist2Post(i) = distance_pts(BB, postPole);
end

% BBs are assigned a position within a cortical row or kinety by using a
% metric that minimizes the distance between each BB and its partner while
% also minimizing the distance between its partner and a plane comprising
% the 3D coordinates of the BB maxima, the anterior pole and the posterior
% pole.

[startPt, endPt] = findLink(x, y, z, antPole, postPole, dist2Ant);
% matrix showing ciliary rows
OriginalBBsMatrix = link(startPt, endPt);

% -------------------------------------------------------------------------

% withLabel contains BBs that belong to a particular ciliary row
% withoutLabel contains BBs that were not assigned to any row
withLabel = intersect(1:num_BBs, OriginalBBsMatrix);
withoutLabel = setdiff(1:num_BBs, withLabel);

% gets rid of ciliary rows that are too short, reassigns the excess BBs to 
% other rows, then resorts these ciliary rows
[~, updatedTraceback] = getIniLabel(OriginalBBsMatrix, withLabel,...
                                    withoutLabel, dist2Ant, x, y, z, ...
                                    antPole, postPole, minBBsInRow);
                                
[f_final, initialMatrix] = reassignBBs(updatedTraceback, x, y, z, ...
                                       antPole, postPole, dist2Ant, ...
                                       dist2Post, withLabel);                      
end
                                     

