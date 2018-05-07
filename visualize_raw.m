% -------------------------------------------------------------------------
% [Ben] 04/30/18 
% Visualizes x-y plane of raw image at depth z. Normalizes the image based
% on the min and max intensity values of the given plane. Have not
% implemented 'cropped' option. 
% ------------------------------------------------------------------------

function visualize_raw(imagepath, z)
channel = 2; % legacy settings
I = readBioImg(imagepath, channel);
layer_z = I(:, :, z);
lo = min(min(layer_z));
hi = max(max(layer_z));
imshow(layer_z, [lo hi]);
end

