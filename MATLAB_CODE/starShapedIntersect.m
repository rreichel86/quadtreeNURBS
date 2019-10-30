function [U] = starShapedIntersect(L1,L2)

%   starShapedIntersect calculates the intersection point of two lines,
%   even if they have to be elongated in order to meet
%   This algorithm is based on the the explanation from Paul Bourke at
%   http://paulbourke.net/geometry/pointlineplane/

% Input:
% L1, L2: the two lines

% Output:
% U: the intersection point

x = [L1(1,1) L2(1,1); L1(1,2) L2(1,2)]; % Sotres the x coordiantes of the lines
y = [L1(2,1) L2(2,1); L1(2,2) L2(2,2)]; % Sotres the y coordiantes of the lines
dx = diff(x);  % Take the differences down each column
dy = diff(y);
den = dx(1)*dy(2)-dy(1)*dx(2);  % Precompute the denominator
ua = (dx(2)*(y(1)-y(3))-dy(2)*(x(1)-x(3)))/den; % Find the "inclination" of each line
ub = (dx(1)*(y(1)-y(3))-dy(1)*(x(1)-x(3)))/den;
xi = x(1)+ua*dx(1); % find x coordinate of intersection
yi = y(1)+ua*dy(1); % find y coordinate of intersection
U = [xi; yi];
end

