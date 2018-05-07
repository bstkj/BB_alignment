% -------------------------------------------------------------------------
% [Ben] 05/01/18 (rewritten by Ben)
% 'first applied a 3D transformation including rotation and translation to
% cell' pg 7. Orients the cell along its antero-posterior axis, so that the
% axis now coincides with the z-axis. Then translates the posterior pole to
% the origin. The BBs coordinates are translated and rotated accordingly.
% -------------------------------------------------------------------------

function [rotated_x, rotated_y, rotated_z] = newCoor(BBs_x, BBs_y, BBs_z, anteriorPole, posteriorPole)

vectorBetweenPoles = anteriorPole - posteriorPole;
r = vrrotvec(vectorBetweenPoles, [0 0 1]); % vector specifying rotation axis and angle
m = vrrotvec2mat(r); % rotation matrix
coordinate_matrix = transpose([BBs_x BBs_y BBs_z]); % each column is an x-y-z coordinate triple
rotated_coordinate_matrix = m*coordinate_matrix;

% anteriorPole and posteriorPole are passed in as 1x3 vectors
% rotated_posteriorPole is a 3x1 vector
rotated_posteriorPole = m*reshape(posteriorPole, 3, 1);

% translate coordinates such that rotated_posteriorPole is shifted to
% origin
rotated_coordinate_matrix = rotated_coordinate_matrix - rotated_posteriorPole;

rotated_x = transpose(rotated_coordinate_matrix(1, :));
rotated_y = transpose(rotated_coordinate_matrix(2, :));
rotated_z = transpose(rotated_coordinate_matrix(3, :));
end

