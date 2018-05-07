% -------------------------------------------------------------------------
% [Ben] 05/02/18
% If loadfrommodel parameter is set to true, then provide modelpath and
% leave imagepath as []. The visualization will be produced based on BB
% position information and CR labels stored in the model specified by the
% modelpath. If loadfrommodel parameter is set to false, then leave
% modelpath as [] and provide imagepath. The BB position information and CR
% labels will be recalculated, saved into a model, and used for the
% visualization.
% -------------------------------------------------------------------------

function visualize_final_CR(loadfrommodel, modelpath, imagepath, slice)

if loadfrommodel
    % should load model struct with the following attributes: cort_x,
    % cort_y, cort_z, oa_x, oa_y, oa_z, antPole, postPole, final_traceback
    data = load(modelpath);
    % these coordinates should already have been rotated and translated as
    % in visualize_initial_CR
    cort_x = data.model.cort_x; oa_x = data.model.oa_x;
    cort_y = data.model.cort_y; oa_y = data.model.oa_y;
    cort_z = data.model.cort_z; oa_z = data.model.oa_z;
    antPole = data.model.antPole; postPole = [0 0 0];
    
    x = vertcat(cort_x, oa_x);
    y = vertcat(cort_y, oa_y);
    z = vertcat(cort_z, oa_z);
    
    % this should reflect the final CR assignments of the cortical BBs
    % (still not including OA BBs because they mess everything up)
    final_traceback = data.model.final_traceback;
else
    channel = 2;
    th = 6;
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
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % This section is for orienting the cell so that the antero-posterior pole
    % is parallel to the z-axis, with the anterior pole lower than the
    % posterior pole on the z-axis.
    [x, y, z] = newCoor(x, y, z, antPole, postPole);
    cort_x = x(1:num_cort); oa_x = x(num_cort+1:end);
    cort_y = y(1:num_cort); oa_y = y(num_cort+1:end);
    cort_z = z(1:num_cort); oa_z = z(num_cort+1:end);
    
    antPole = [0 0 norm(antPole - postPole)];
    postPole = [0 0 0];
    
    x = vertcat(cort_x, oa_x);
    y = vertcat(cort_y, oa_y);
    z = vertcat(cort_z, oa_z);
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % This section is for saving the trained model to the models folder. 
    % Still need to save hyperparameters such as minBBsInRow, th, etc.
    model.cort_x = cort_x; model.oa_x = oa_x;
    model.cort_y = cort_y; model.oa_y = oa_y;
    model.cort_z = cort_z; model.oa_z = oa_z;
    model.final_traceback = final_traceback;
    [~, name, ~] = fileparts(imagepath);
    model.name = name;
    save(sprintf('models/%s.mat', name), 'model');
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

num_rows = size(final_traceback, 1);
if strcmp(slice, 'XYZ')
    K = convhull(x, y, z);
    % plot mesh of convex hull of BB positions
    trimesh(K, x, y, z, 'FaceColor', 'red', 'FaceAlpha', 0.1, ...
        'EdgeColor', 'red', 'EdgeAlpha', 0);
    hold on
    % plot cortical BBs as red filled circles
    scatter3(cort_x, cort_y, cort_z, 80, 'r.');
    % plot oral apparatus BBs as green filled circles
    scatter3(oa_x, oa_y, oa_z, 80, 'g.');
    % plot black straight line between the anterior and posterior poles
    plot3([postPole(1) antPole(1)], [postPole(2) antPole(2)], ...
        [postPole(3) antPole(3)], '-', 'Color', 'black', 'LineWidth', 3);
    
    % color consecutive ciliary rows with different colors for contrast
    colors = {'yellow', 'magenta', 'cyan', 'blue'};
    xy_centroids = zeros(num_rows, 1);
    for i = 1:num_rows
        row = final_traceback(i, :);
        row = row(row ~= 0);
        xy_centroids(i) = atan2(mean(cort_y(row)), mean(cort_x(row)));
    end
    [~, color_order] = sort(xy_centroids, 'ascend');
    for i = 1:num_rows
        row = final_traceback(color_order(i), :);
        row = row(row ~= 0);
        color = colors{mod(i, 4) + 1};
        plot3(cort_x(row), cort_y(row), cort_z(row), 'Color', color, ...
            'LineStyle', '-', 'LineWidth', 1);
    end
    daspect([1 1 (0.125/0.3)]);
    zlabel(slice(3));
elseif strcmp(slice, 'XY')
    K = convhull(x, y);
    fill(x(K), y(K), [1 1 1]);
    hold on
    scatter(cort_x, cort_y, 80, 'r.'); 
    for i = 1:num_rows
        row = final_traceback(i, :);
        row = row(row ~= 0);
        plot(cort_x(row), cort_y(row), 'b-');
    end
elseif strcmp(slice, 'XZ')
    K = convhull(x, z);
    fill(x(K), z(K), [1 1 1]);
    hold on
    scatter(cort_x, cort_z, 80, 'r.');
    for i = 1:num_rows
        row = final_traceback(i, :);
        row = row(row ~= 0);
        plot(cort_x(row), cort_z(row), 'b-');
    end
    daspect([1 (0.125/0.3) 1]);
elseif strcmp(slice, 'YZ')
    K = convhull(y, z);
    fill(y(K), z(K), [1 1 1]);
    hold on
    scatter(cort_y, cort_z, 80, 'r.');
    for i = 1:num_rows
        row = final_traceback(i, :);
        row = row(row ~= 0);
        plot(cort_y(row), cort_z(row), 'b-');
    end
    daspect([1 (0.125/0.3) 1]);
end
title(slice);
xlabel(slice(1));
ylabel(slice(2));
hold off
end

