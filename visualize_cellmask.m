% -------------------------------------------------------------------------
% [Ben] 01/02/18 (written by Ben)
% Visualize cell mask for given .nd2 file in /nd2_files. This function
% takes in a path to a .nd2 file, and takes about 8s to run. Produces
% figure 2.2 b. 
% -------------------------------------------------------------------------

function visualize_cellmask(imagepath)
channel = 2; % legacy settings
I = readBioImg(imagepath, channel);
[cell, ~] = identifyCell(I);
imshow(cell);
end

