function [ori] = starShapedCheckOrientation(vertices)

%   starShapedCheckOrientation checks the orientation of an polygon
%   Returns -1 if polygon is in clockwise direction and 1 otherwise

% Input:
% vertices
% Output:
% ori: -1 = clockwise, 1 = counter-clockwise, 0 = all vertices on a line

verticesFirsttoEnd = [vertices,vertices(:,1)]; % copies the first point to the end of the matrix

% Find Orientation
if round(polyarea(vertices(1,:), vertices(2,:)), 5) == 0
    ori = 0; % all points on a line
else
    % Apply the sum over edges rule: Sum (x2-x1)*(y2-y1) for all edges, of
    % the sum is positive the polyogn is given in clockwise direction,
    % if the sum is negative, in counter-clockwise. If the sum is 1, it
    % means that all th given vertices lie on a line
    sum = 0;
    for p = 1:length(vertices)
        x1 = verticesFirsttoEnd(1,p);
        x2 = verticesFirsttoEnd(1,p+1);
        y1 = verticesFirsttoEnd(2,p);
        y2 = verticesFirsttoEnd(2,p+1);
        sum = sum + (x2-x1)*(y2+y1);
    end
    if sum < 0
        ori = 1;
    else
        ori = -1;
    end
end

