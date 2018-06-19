% -------------------------------------------------------------------------
% [Ben] 05/28/18
% Basically a modified 'main.m'. modelname is the name (excluding the .mat
% extension) of the saved model file. If modelname is set to [], then the
% default value used will be the name of the .nd2 file (excluding the .nd2
% extension). In addition to saving model parameters to a .mat file, now
% also saves model parameters to a .xlsx file.
% Should separate the parameter saving from the parameter extraction.
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
    modelname = name;
end
model.name = modelname;

% for saving model parameters to a .mat file in models/
try
    save(sprintf('models/%s.mat', model.name), 'model');
catch exception
    msgText = getReport(exception);
    warning(msgText);
    pause
end

% for saving model parameters to an .xlsx file in models_xlsx/
try
    T_cort = table(cort_x, cort_y, cort_z, ...
        'VariableNames', {'cort_x' 'cort_y' 'cort_z'});
    T_oa = table(oa_x, oa_y, oa_z, ...
        'VariableNames', {'oa_x' 'oa_y' 'oa_z'});
    T_antPole = table(antPole(1), antPole(2), antPole(3), ...
        'VariableNames', {'antPole_x' 'antPole_y' 'antPole_z'});
    % note that transpose is saved to .xlsx file, but for .mat file the
    % matrix saved is not transposed
    T_finalTraceback = array2table(transpose(final_traceback), ...
        'VariableNames', cellfun(@(x) sprintf('R%d', x), ...
        num2cell(1:size(final_traceback, 1)), 'UniformOutput', false));
    filename = sprintf('models_xlsx/%s.xlsx', model.name);
    % this information can later be extracted from the .xlsx file by
    % readtable(filename, 'Sheet', <sheet_number>)
    writetable(T_cort, filename, 'Sheet', 1);
    writetable(T_oa, filename, 'Sheet', 2);
    writetable(T_antPole, filename, 'Sheet', 3);
    writetable(T_finalTraceback, filename, 'Sheet', 4);
catch exception
    msgText = getReport(exception);
    warning(msgText);
    pause
end
end