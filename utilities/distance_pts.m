% -------------------------------------------------------------------------
% [Ben] 12/07/17
% Calculates and returns the Euclidean distance between 2 points.
% -------------------------------------------------------------------------

function d=distance_pts(point1, point2)
% unit of d is micron

% x1=point1(1)*0.125;
% y1=point1(2)*0.125;
% z1=point1(3)*0.3;
% 
% x2=point2(1)*0.125;
% y2=point2(2)*0.125;
% z2=point2(3)*0.3;

x1=point1(1)*0.125;
y1=point1(2)*0.125;
z1=point1(3)*0.125;

x2=point2(1)*0.125;
y2=point2(2)*0.125;
z2=point2(3)*0.125;

d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);

end