% -------------------------------------------------------------------------
% [Ben] 05/03/18
% Executes visualization tools that do not take that long to run, and so do
% not have the option to load variables directly from a saved model.
% -------------------------------------------------------------------------

function demo001()

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

% xy slice of raw image at given depth, z
z = 5;
visualize_raw(imagepath, z);
pause

% cellmask for cropping raw image
visualize_cellmask(imagepath);
pause

% maximum intensity projection of raw image
visualize_MIP(imagepath);
pause

% identified BBs and convex hull built from their positions
slice = 'XYZ';
visualize_BB(imagepath, slice);
pause

% identified BBs, and initial ciliary row assignments for cortical BBs (BBs
% that are drawn as red filled circles, as opposed to green filled circles
% which denote oral apparatus BBs)
visualize_initial_CR(imagepath, slice);
pause

end

