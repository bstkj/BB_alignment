% -------------------------------------------------------------------------
% [Ben] 04/30/18 (written by Ben)
% Visualize maximum intensity projection for given .nd2 file in /nd2_files.
% Takes about 8s to run. Produces figure 2.2 a. 
% -------------------------------------------------------------------------

function visualize_MIP(imagepath)
channel = 2; % legacy settings
I = readBioImg(imagepath, channel);
MIP = max(I, [], 3); % return maxes along dim 3 (depth)
lo = min(min(MIP));
hi = max(max(MIP));
imshow(MIP, [lo hi]);
end