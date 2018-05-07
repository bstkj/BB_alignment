% -------------------------------------------------------------------------
% [Ben] 12/08/17
% Takes in the valid anterior-posterior BB pairs and returns a matrix
% showing the reconstructed ciliary rows. Each row of the matrix contains a
% sequence of BBs (padded with 0's at the end) that make up a ciliary row.
% -------------------------------------------------------------------------

function [traceback] = link(startPt, endPt)

onlyInStart = [];
% find the BBs that are not anterior partners
for i = 1:length(startPt)
    if isempty(find(endPt == startPt(i), 1))
        onlyInStart = [onlyInStart, startPt(i)];
    end
end
onlyInEnd = [];
% find the BBs that are not posterior partners
for i = 1:length(endPt)
    if isempty(find(startPt == endPt(i), 1))
        onlyInEnd = [onlyInEnd, endPt(i)];
    end
end

% traceback(:, 1) = onlyInStart; % traceback is the column vector of onlyInStart
traceback = transpose(onlyInStart);

for i = 1:length(onlyInStart)
    j = 1;
    % essentially building the ciliary rows
    while isempty(find(onlyInEnd == traceback(i, j), 1))
        idx = find(startPt == traceback(i, j));
        traceback(i, j+1) = endPt(idx);
        j = j+1;
    end
end

end