% -------------------------------------------------------------------------
% [Ben] 05/03/18
% Shows how to use model2report tool.
% -------------------------------------------------------------------------

function demo004()

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

if exist('models/003.mat', 'file')
    model2report('models/003.mat')
else
    disp('Model file "003.mat" does not exist, so we are training it now.');
    imagepath = 'nd2_files/003.nd2';
    modelname = [];
    train_model(imagepath, modelname);
    disp('Model file "003.mat" has been produced.');
    model2report('models/003.mat')
end
end


