% -------------------------------------------------------------------------
% [Ben] 1/30/18
% For each ciliary row, calculates the average BB position, then sums the
% distances between each BB in the ciliary row and the plane defined by the
% average BB position, the anterior pole, and the posterior pole. Totals
% these sums across all ciliary rows and returns the result.
% -------------------------------------------------------------------------

function f3 = constraint3(mtx, antPole, postPole, x, y, z)
numRows = size(mtx, 1);
row_counts = sum(mtx ~= 0, 2);
f3 = 0;
for row = 1:numRows
    numCols = row_counts(row);
    currRow = mtx(row, 1:numCols);
    pt1 = [mean(x(currRow)), mean(y(currRow)), mean(z(currRow))];
    for col = 1:numCols
        bb = mtx(row, col);
        f3 = f3 + distance_pt2plane([x(bb), y(bb), z(bb)], pt1, antPole, postPole); 
    end
end
end
