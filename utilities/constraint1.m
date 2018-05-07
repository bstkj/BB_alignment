% -------------------------------------------------------------------------
% [Ben] 1/30/18
% Takes in a traceback matrix (matrix whose rows correspond to the ciliary
% rows of a cell, with each element in each row representing a particular
% basal body) and returns the product of the sample variance in the no. of
% BB's per ciliary row, multiplied by the no. of ciliary rows.
% -------------------------------------------------------------------------

function f1 = constraint1(mtx)
row_counts = sum(mtx ~= 0, 2);
f1 = sum((row_counts - mean(row_counts)).^2);
end