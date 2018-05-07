% -------------------------------------------------------------------------
% [Ben] 1/30/18
% Takes in traceback matrix as well as x, y, z coordinates of all detected
% BBs. For each ciliary row, calculates the max distance to be found
% between two adjacent member BBs. Returns the sum of these 'max distances'
% across all ciliary rows. 
% -------------------------------------------------------------------------

function f2 = constraint2(mtx, x, y, z)
x = x*0.125;
y = y*0.125;
z = z*0.3;
numRows = size(mtx, 1);
row_counts = sum(mtx ~= 0, 2);
f2 = 0;
for row = 1:numRows
    max_d = 0;
    % mtx will not have rows whose number of non-zero entries < 2
    for col = 2:row_counts(row)
        i = mtx(row, col);
        j = mtx(row, col-1); 
        d = sum(([x(i), y(i), z(i)] - [x(j), y(j), z(j)]).^2);
        if d > max_d
            max_d = d;
        end
    end
    f2 = f2 + sqrt(max_d);
end
end