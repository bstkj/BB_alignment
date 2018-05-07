% -------------------------------------------------------------------------
% [Ben] 12/08/17
% Returns the value of the objective function we are trying to minimize in
% our BB reassignment algorithm. 
% -------------------------------------------------------------------------

function f=constraints(Mtx, antPole, postPole, dist2Ant, dist2Post, x, y, z, u)
f1 = constraint1(Mtx);
f2 = constraint2(Mtx, x, y, z);
f3 = constraint3(Mtx, antPole, postPole, x, y, z);
f4 = constraint4(Mtx, dist2Ant, dist2Post);
f5 = constraint5(Mtx, x, y, z);
f = dot([f1, f2, f3, f4, f5], u);
end