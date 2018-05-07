% -------------------------------------------------------------------------
% [Ben] 05/03/18
% Executes visualization tools that take that long to run if run from
% scratch (because the model of the ciliary rows needs to be trained
% first). Also shows how to first train and save a model so that it can be
% used by various visualization tools.
% -------------------------------------------------------------------------

function demo002()

[~, name, ~] = fileparts(pwd);
if ~strcmp(name, 'BB_alignment')
    cd('../'); % assuming you were in some subfolder of BB_alignment
    [~, name, ~] = fileparts(pwd);
    if ~strcmp(name, 'BB_alignment') 
        % if you are still not in BB_alignment after moving one level up,
        % then you are in the wrong place
        disp('Please run this function within the BB_alignment folder');
        return
    end
    setup;
end

% Train a model on one of the .nd2 images in the folder nd2_files, and save
% it
% THIS TAKES SOME TIME TO RUN (On my slow laptop it takes about 50 minutes)
imagepath = 'nd2_files/003.nd2';
modelname = [];
% because we set modelname to [], by default our model is going to be saved
% as 003.mat in the models folder
train_model(imagepath, modelname);
disp('Trained model is saved as 003.mat in the models folder.');
pause;

% Visualize the final ciliary row assignments that were trained.
loadfrommodel = true;
modelpath = 'models/003.mat'; 
% imagepath is set to [] because we are taking information directly from
% the trained model and not from the raw image
imagepath = [];
% slice is set to 'XYZ' because we want to visualize the ciliary row
% assignments in 3D
slice = 'XYZ';
visualize_final_CR(loadfrommodel, modelpath, imagepath, slice)
pause;
end

