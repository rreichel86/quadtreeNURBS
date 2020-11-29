function [Quadtree] = Star_Shape(Quadtree,NURBS,leaves)
% Star_Shape: Check if Quadtree leaves are star-shaped. Only the Quadtree 
% leaves that are splitted by the NURBS curve are checked. 
%
% INPUT: 
% Quadtree ------------------- Quadtree data structure
% NURBS definition 
% NURBS.degree --------------------- NURBS degree
% NURBS.knots ---------------------- NURBS knot vector
% NURBS.controlPoints -------------- NURBS control points
% NURBS.weights -------------------- NURBS weights
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
        [star_shaped, Quadtree] = starAlgorithm1(Quadtree,leaves(i));
        if star_shaped == 1
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
    [Quadtree] = Decompose_helper(Quadtree,NURBS,leaves(i));
    % get children pointers of current Quad 
    newLeaves(4*i-3:4*i) = Quadtree.Node{leaves(i),1}{11,1};
end

% chech if new Quadtree leaves are shar-shaped
[Quadtree] = Star_Shape(Quadtree,NURBS, newLeaves);


