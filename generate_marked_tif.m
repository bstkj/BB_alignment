% -------------------------------------------------------------------------
% [Ben] 05/02/18
% Takes in path to a raw .nd2 image, picks out the BBs in each plane, then
% generates a tif marking the cortical BBs and the oral apparatus BBs
% detected in each xy slice of the .nd2 image. 
% -------------------------------------------------------------------------

function generate_marked_tif(imagepath, tifname)
channel = 2;
th = 6;
I = readBioImg(imagepath, channel);

masks = create_BB_masks(I, th);

[rows, cols, layers] = size(I);
output = zeros(rows, cols, 3, layers); % 3 channels

% for marking cortical BBs in current plane
redpatch = cat(3, [1 1 1; 1 1 1; 1 1 1], [0 0 0; 0 0 0; 0 0 0], ...
    [0 0 0; 0 0 0; 0 0 0]);
% for marking oral BBs in current plane
greenpatch = cat(3, [0 0 0; 0 0 0; 0 0 0], [1 1 1; 1 1 1; 1 1 1], ...
    [0 0 0; 0 0 0; 0 0 0;]);
% for marking cortical/oral BBs in adjacent planes to current plane
bluepatch = cat(3, [0 0 0; 0 0 0; 0 0 0], [0 0 0; 0 0 0; 0 0 0], ...
    [1 1 1; 1 1 1; 1 1 1]);

% marks the BBs on the raw image matrix
for z = 1:layers
    I_layer = I(:, :, z);
    mx = max(max(I_layer));
    mn = min(min(I_layer));
    norm_I_layer = (I_layer - mn)/(mx - mn);
    [cort_y, cort_x] = find(masks(:, :, z) == 1);
    [oa_y, oa_x] = find(masks(:, :, z) == 2);
    
        
    for c = 1:3
        output(:, :, c, z) = norm_I_layer;
    end
    
    yd = [];
    xd = []; 
    if z-1 > 0
        [yd, xd] = find(masks(:, :, z-1) > 0);
        if z-2 > 0
            [yd2, xd2] = find(masks(:, :, z-2) > 0);
            yd = vertcat(yd, yd2);
            xd = vertcat(xd, xd2);
        end
    end
    
    yu = [];
    xu = [];
    if z+1 <= layers
        [yu, xu] = find(masks(:, :, z+1) > 0);
        if z+2 <= layers
            [yu2, xu2] = find(masks(:, :, z+2) > 0);
            yu = vertcat(yu, yu2);
            xu = vertcat(xu, xu2);
        end
    end
    
    for i = 1:length(yu)
        output((yu(i)-1):(yu(i)+1), (xu(i)-1):(xu(i)+1), :, z) = bluepatch;
    end
    
    for i = 1:length(yd)
        output((yd(i)-1):(yd(i)+1), (xd(i)-1):(xd(i)+1), :, z) = bluepatch;
    end
    
    for i = 1:length(cort_y)
        output((cort_y(i)-1):(cort_y(i)+1), (cort_x(i)-1):(cort_x(i)+1), :, z) = redpatch;
    end
    
    for i = 1:length(oa_y)
        output((oa_y(i)-1):(oa_y(i)+1), (oa_x(i)-1):(oa_x(i)+1), :, z) = greenpatch;
    end
end

if isempty(tifname)
    [~, name, ~] = fileparts(imagepath);
else
    name = tifname;
end
target_filename = ['generated_tifs' filesep name '.tif'];

% saves the marked image matrix as a .tif file in the generated_tifs folder
for z = 1:layers
    try
        imwrite(output(:, :, :, z), target_filename, 'WriteMode', 'append');
    catch
        imwrite(output(:, :, :, z), target_filename);
    end
end
end

