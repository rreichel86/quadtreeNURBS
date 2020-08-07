function [Quadtree]=Star_Shape(Quadtree, controlPoints, knots, weights, degree, Boundary, leaves)
% Star_Shape: Check if Quadtree leaves are star-shaped. Only the Quadtree 
% leaves that are splitted by the NURBS curve are checked. 
%
% INPUT: 
% Quadtree ------------------- Quadtree data structure
% controlPoints -------------- NURBS control points 
% knots ---------------------- NURBS knot vector
% weights -------------------- NURBS weights
% degree --------------------- NURBS degree
% Boundary
% leaves --------------------- filtered Quadtree leaves (optional)
%
% OUTPUT: 
% Quadtree ------------------- updated Quadtree data structure
%
% -------------------------------------------------------------------------

% Find Quadtree leaves
if ~exist('leaves','var')
leaves = Quadtree.findleaves();
end 

% Filter Quadtree leaves
for i = 1:length(leaves)
    intersections=Quadtree.Node{leaves(i),1}{5,1};
    % Quadtree leaves with 0 or 1 intersection
    if isempty(intersections) || length(intersections) == 1
        leaves(i) = 0;      
    else  % Quadtree leaves with 2 intersections 
        % use algorithm 1
        % check if the 2 Quad subregions are star-shaped
        [Quadtree] = starAlgorithm1(Quadtree,leaves(i));
        if ~isempty(Quadtree.Node{leaves(i),1}{12,1}) && ~isempty(Quadtree.Node{leaves(i),1}{13,1})
            leaves(i) = 0; % delete Quad from the list
            continue
        end
        % if algorithm 1 fails use algorithm 2
        % check if the 2 Quad subregions are star-shaped
        [Quadtree] = starAlgorithm2(Quadtree,leaves(i));
        if ~isempty(Quadtree.Node{leaves(i),1}{12,1}) && ~isempty(Quadtree.Node{leaves(i),1}{13,1})
            leaves(i) = 0; % delete Quad from the list
            continue
        end 
    end
end
% update leaves 
leaves = leaves(leaves~=0);

if isempty(leaves)
    return
end     

% subdivide remainig leaves
newLeaves = zeros(1, 4*length(leaves));
for i = 1:length(leaves)
    [Quadtree] = Decompose_Star(Quadtree, controlPoints, ...
    knots, weights, degree,leaves(i),Boundary);
    % get children pointers of current Quad 
    newLeaves(4*i-3:4*i) = Quadtree.Node{leaves(i),1}{11,1};
end

% chech if new Quadtree leaves are shar-shaped
[Quadtree] = Star_Shape(Quadtree,controlPoints, knots, weights, degree, Boundary, newLeaves);


