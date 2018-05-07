% -------------------------------------------------------------------------
% [Ben] 1/31/18
% Moves pt from its original location in oldMtx to row newRow, then sorts
% the row.
% -------------------------------------------------------------------------

function newMtx = moveBB(oldMtx, pt, newRow, dist2Ant)
% we assume the following invariants on the input arguments:
% 1) oldMtx does not have any trailing columns of all zeros on the right
% 2) the BB indices in each row of oldMtx are already sorted in order of
% decreasing distance from the anterior pole

[oldRow, oldCol] = find(oldMtx == pt);

if newRow == oldRow
    % inserting BB into the row that it is already assigned will have no
    % change
    newMtx = oldMtx;
else
    row_counts = sum(oldMtx ~= 0, 2);
    numBBs_oldRow = row_counts(oldRow); % no. of BBs in row 'oldRow'
    numBBs_newRow = row_counts(newRow); % no. of BBs in row 'newRow'

    if oldMtx(newRow, end) == 0
        % newMtx keeps the same dimensions as oldMtx
        newMtx = oldMtx;
    else
        % newMtx has one more column than oldMtx
        newMtx = [oldMtx, zeros(size(oldMtx, 1))];
    end
    
    % insert pt into row 'newRow'
    bbs = [oldMtx(newRow, 1:numBBs_newRow), pt];
    [~, I] = sort(dist2Ant(bbs), 'descend');
    newMtx(newRow, 1:(numBBs_newRow + 1)) = bbs(I);

    % remove pt from row 'oldRow'
    if oldCol ~= numBBs_oldRow
        newMtx(oldRow, oldCol:(numBBs_oldRow - 1)) = oldMtx(oldRow, (oldCol+1):numBBs_oldRow);
    end
    newMtx(oldRow, numBBs_oldRow) = 0;
end