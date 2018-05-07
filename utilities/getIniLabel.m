% -------------------------------------------------------------------------
% [Ben] 2/1/18
% Returns an updated traceback matrix, together with the corresponding
% label vector. The traceback matrix consists of sequences of BBs (padded
% with trailing 0's) in each row, tracing out a ciliary row. The label
% vector indicates which BBs belong to the same ciliary row. The updated
% traceback matrix may have fewer rows than the original traceback matrix.
% In effect, this function deletes ciliary rows with too few BBs and
% reassigns the BBs from deleted ciliary rows to other ciliary rows.
% -------------------------------------------------------------------------


function [iniLabel, updatedTraceback] = getIniLabel(mtx, withLabel, withoutLabel, dist2Ant, x, y, z, antPole, postPole, minBBsInRow)
% we assume the following invariants on the input arguments:
% 1) the BB indices in each row of mtx are already sorted in order of
% decreasing distance from the anterior pole.

[numRows, ~] = size(mtx);
CRs = 1:numRows;
row_counts = sum(mtx ~= 0, 2);
legalRows = CRs(row_counts >= minBBsInRow);
legal_row_counts = row_counts(legalRows);
numLegalRows = length(legalRows);

% numRows is the no. of ciliary rows detected by findLink.m and link.m
numBBs = length(withLabel) + length(withoutLabel);
label = zeros(numBBs, 1);

for i = 1:numLegalRows
    row = legalRows(i);
    bbs = mtx(row, 1:row_counts(row));
    label(bbs) = i;
end

if numLegalRows ~= numRows
    % update traceback matrix by only including legal rows
    legalRows_mtx = mtx(legalRows, :);
    nonLegal_BBs = setdiff(find(label == 0), withoutLabel);
    num_nonLegal_BBs = length(nonLegal_BBs);
    for i = 1:num_nonLegal_BBs
        bb = nonLegal_BBs(i);
        % assign BBs not in legal rows to a legal row based on their
        % distance to the 'average planes' corresponding to each legal row
        [~, rowIdx] = dist2plane(legalRows_mtx, x, y, z, bb, antPole, postPole);
        % legal_row_counts changes every iteration
        legal_row_counts(rowIdx) = legal_row_counts(rowIdx) + 1;
        % legalRows_mtx changes size every iteration
        legalRows_mtx(rowIdx, legal_row_counts(rowIdx)) = bb;
        % update 'label' vector
        label(bb) = rowIdx;
    end
    updatedTraceback = sort_d2ant(legalRows_mtx, dist2Ant);
else
    % no need to sort in this case
    updatedTraceback = mtx;
end
iniLabel = label;
end