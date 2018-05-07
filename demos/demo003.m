% -------------------------------------------------------------------------
% [Ben] 05/03/18
% Shows how to use generated_marked_tif tool.
% -------------------------------------------------------------------------

function demo003()

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

imagepath = 'nd2_files/003.nd2';
% because we set tifname to [], by default our model is going to be saved
% as 003.tif in the models folder
tifname = [];
generate_marked_tif(imagepath, tifname)
end

