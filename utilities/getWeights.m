% -------------------------------------------------------------------------
% [Ben] 12/08/17
% Calculates the normalizing weights for use in the objective function that
% we will try to minimize as we try out various BB reassignments. Returns
% the weights as a vector.
% -------------------------------------------------------------------------


function u=getWeights(oriMtx, antPole, postPole, x, y, z, dist2Ant, dist2Post)

% pg. 20 of Jingyi's paper. 
f1_ini = constraint1(oriMtx);
u1_ini=1/f1_ini;
f2_ini = constraint2(oriMtx, x, y, z);
u2_ini = 1/f2_ini;
f3_ini = constraint3(oriMtx, antPole, postPole, x, y, z);
u3_ini = 1/f3_ini;
f4_ini = constraint4(oriMtx, dist2Ant, dist2Post);
u4_ini = 1/f4_ini;
f5_ini = constraint5(oriMtx, x, y, z);
u5_ini = 1/f5_ini;
% rescale u
% what is the purpose of this?
u=u1_ini+u2_ini+u3_ini+u4_ini+u5_ini;
scale=1/u;
u(1)=u1_ini*scale;
u(2)=u2_ini*scale;
u(3)=u3_ini*scale;
u(4)=u4_ini*scale;
u(5)=u5_ini*scale;
end