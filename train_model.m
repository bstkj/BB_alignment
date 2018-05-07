% -------------------------------------------------------------------------
% [Ben] 05/02/18
% Basically a modified 'main.m'. modelname is the name (excluding the .mat
% extension) of the saved model file. If modelname is set to [], then the
% default value used will be the name of the .nd2 file (excluding the .nd2
% extension).
% -------------------------------------------------------------------------

function train_model(imagepath, modelname)
channel = 2; % legacy settings
th = 6; % legacy settings
I = readBioImg(imagepath, channel);

[cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, ~, ~, ~, ~] = getBBIdx(I, th);
x = vertcat(cort_x, oa_x);
y = vertcat(cort_y, oa_y);
z = vertcat(cort_z, oa_z);

OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);

[pole1, pole2] = findPoles(cort_x, cort_y, cort_z);
d1 = distance_pts(pole1, OA);
d2 = distance_pts(pole2, OA);
% anterior pole is closer to OA region
if d1 < d2
    antPole = pole1;
    postPole = pole2;
else
    antPole = pole2;
    postPole = pole1;
end

minBBsInRow = 10; % legacy settings
% repeated calculation of antPole and postPole within ... remove
% redundancy later
% note that we are going to calculate final_traceback based on
% cortical BBs and are excluded OA BBs for now.
[~, final_traceback] = getCiliaryRows(cort_x, cort_y, cort_z, ...
    antPole, postPole, minBBsInRow);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Rotation and translation of BB coordinates
[x, y, z] = newCoor(x, y, z, antPole, postPole);
num_cort = length(cort_x);
cort_x = x(1:num_cort);
cort_y = y(1:num_cort);
cort_z = z(1:num_cort);
oa_x = x(num_cort+1:end);
oa_y = y(num_cort+1:end);
oa_z = z(num_cort+1:end);

antPole = [0 0 norm(antPole - postPole)];
% postPole = [0 0 0];
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

model.cort_x = cort_x; model.oa_x = oa_x;
model.cort_y = cort_y; model.oa_y = oa_y;
model.cort_z = cort_z; model.oa_z = oa_z;
% do not need to store postPole, since it is always at the origin after the
% rotation and translation process
model.antPole = antPole;
model.final_traceback = final_traceback;

% should eventually turn to using try-catch statements for error handling
if isempty(modelname)
    [~, name, ~] = fileparts(imagepath);
    model.name = name;
else
    model.name = modelname;
end

save(sprintf('models/%s.mat', model.name), 'model');
end

