% -------------------------------------------------------------------------
% [Ben] 1/31/18
% For each ciliary row, calculates the minimum distance among all member
% BBs to the anterior pole, and the minimum distance to the posterior pole.
% Takes the sum of these 2 distances, adds them up across all ciliary rows
% and returns the total. 
% -------------------------------------------------------------------------

function f4 = constraint4(mtx, dist2Ant, dist2Post)
numRows = size(mtx, 1);
row_counts = sum(mtx ~= 0, 2);
f4 = 0;
for row = 1:numRows
    numCols = row_counts(row);
    currRow = mtx(row, 1:numCols);
    min_ant_d = Inf;
    min_post_d = Inf;
    for col = 1:numCols
        bb = currRow(col);
        ant_d = dist2Ant(bb);
        if ant_d < min_ant_d
            min_ant_d = ant_d;
        end
        post_d = dist2Post(bb);
        if post_d < min_post_d
            min_post_d = post_d;
        end
    end
    f4 = f4 + min_ant_d + min_post_d;
end
end