% -------------------------------------------------------------------------
% [Ben] 05/28/18
% Run batch jobs for train_model in parallel. imagepaths / modelnames
% should be a cell-array of character vectors /strings. If modelnames is
% set to {}, then the default value for modelnames will be used.
% Have to do more research on parallelism in matlab.
% -------------------------------------------------------------------------

function train_model_batch(imagepaths, modelnames)

n = length(imagepaths);

if isempty(modelnames)
    % modelnames is initialized to nx1 cell array of empty matrices
    modelnames = cell(n, 1);
end

try 
    parfor i = 1:n
        imagepath = imagepaths{i};
        modelname = modelnames{i};
        train_model(imagepath, modelname);
    end
catch exception
    msgText = getReport(exception);
    warning(msgText);
    pause
end

% get current parallel pool and shut it down
poolobj = gcp('nocreate');
delete(poolobj);

end

