% -------------------------------------------------------------------------
% [Ben] 1/31/18
% For each ciliary row, calculates the vectors between adjacent BB's, 
% computes the angles between consecutive pairs of vectors, and sums up
% these angular differences. Then adds these sums across all ciliary rows
% and returns the total.
% -------------------------------------------------------------------------

function f5 = constraint5(mtx, x, y, z)
x = x*0.125;
y = y*0.125;
z = z*0.3;

numRows=size(mtx, 1);
% angle in radians
f5 = 0;
for row=1:numRows
    currRow=mtx(row, :);
    currRow(currRow == 0) = [];
    BBCount=length(currRow);
    for bb = 2:(BBCount-1)
        v1=[x(currRow(bb)) - x(currRow(bb-1)), y(currRow(bb)) - y(currRow(bb-1)), ...
            z(currRow(bb)) - z(currRow(bb-1))];
        v2=[x(currRow(bb + 1)) - x(currRow(bb)), y(currRow(bb + 1)) - y(currRow(bb)), ...
            z(currRow(bb + 1)) - z(currRow(bb))];
        angle = atan2(norm(cross(v1, v2)), dot(v1, v2));
        f5 = f5+angle;
    end
end
end