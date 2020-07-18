function [numel] = countElements(Quadtree)
% countElements: determine total number of elements

% get leaves
l = Quadtree.findleaves();

% loop over the leaves
% check leaf for intersection point
% The leaf is splitted in 2 elements if it has 2 intersection points
numel = 0;
for i = 1:length(l)
    intersections = Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        numel = numel + 1;
    else
        numel = numel + 2;
    end
end
end
