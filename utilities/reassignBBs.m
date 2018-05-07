% -------------------------------------------------------------------------
% [Ben] 2/2/18 (added by Ben)
% Performs the greedy reassignment of BBs to most appropriate ciliary rows.
% Works on an updatedTraceback matrix that already only contains legal
% rows (the output of getIniLabel). Returns the F-score of the final
% ciliary row configuration arrived at, as well as the corresponding
% traceback matrix.
% -------------------------------------------------------------------------


function [f_final, initialMatrix] = reassignBBs(updatedTraceback, x, y, z, antPole, postPole, dist2Ant, dist2Post, withLabel)

coeffs = getWeights(updatedTraceback, antPole, postPole, x, y, z, dist2Ant, dist2Post);
f_old = constraints(updatedTraceback, antPole, postPole, dist2Ant, dist2Post, x, y, z, coeffs);

initialMatrix = updatedTraceback;
[numRow, ~] = size(initialMatrix);
numBB_withLabel = length(withLabel);

% pg. 20 of Jingyi's paper. Mtx is a matrix with n rows and m columns,
% where n is the no. of valid ciliary rows in the cell, and m is the no. of
% BBs that have been assigned to a ciliary row. Mtx(a, b) saves the
% difference between the F-score of the current BB assignment and the
% assignment if we assigned BB b to row a.
minMtxVal = -Inf;
while minMtxVal < 0
    Mtx = zeros(numRow, numBB_withLabel);
    for row = 1:numRow
        for bb = 1:numBB_withLabel
            newMtx = moveBB(initialMatrix, withLabel(bb), row, dist2Ant);
            f = constraints(newMtx, antPole, postPole, dist2Ant, dist2Post, ...
                            x, y, z, coeffs);
            Mtx(row, bb) = f - f_old;
        end
    end
    
    % do not suppress output so we can observe the convergence 
    minMtxVal = min(min(Mtx))
    f_old = f_old + minMtxVal;
    
    % if there are several entries that all have the same min value, then
    % just pick the first pair of indices
    [minValRowIdx, minValBBIdx] = find(Mtx == minMtxVal, 1);
    initialMatrix = moveBB(initialMatrix, withLabel(minValBBIdx), minValRowIdx, dist2Ant);
end

f_final = f_old;
end

