% -------------------------------------------------------------------------
% [Ben] 1/31/18
%  Sorts BBs in the same ciliary row according to their distances to the
%  anterior pole. Returns the resultant sorted traceback matrix (with only
%  legal rows).
% -------------------------------------------------------------------------

function sortedMtx = sort_d2ant(oriMtx, dist2Ant)
% we assume the following invariants on the input arguments:
% 1) oriMtx does not have any trailing columns of all zeros on the right

[numRows, numCols] = size(oriMtx);
row_counts = sum(oriMtx ~= 0, 2);
sortedMtx = zeros(numRows, numCols);
for row = 1:numRows
    numBBs_row = row_counts(row);
    bbs = oriMtx(row, 1:numBBs_row);
    [~, I] = sort(dist2Ant(bbs), 'descend');
    sortedMtx(row, 1:numBBs_row) = bbs(I);
end
end