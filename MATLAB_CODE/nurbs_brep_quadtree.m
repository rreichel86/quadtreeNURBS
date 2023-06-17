function [Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary,dbg)
% nurbs_brep_quadree: quadtree decomposition 
% 
% INPUT: 
% k_min ---------------------- minimum level of decomposition
% NURBS definition
% NURBS.degree --------------- NURBS degree
% NURBS.knots ---------------- NURBS knot vector
% NURBS.controlPoints -------- NURBS control points
% NURBS.weights -------------- NURBS weights
% Boundary ------------------- Outer boundary 
% 
% OUTPUT:
% Quadtree ------------------- Quadtree data structure 
% 
% -------------------------------------------------------------------------

if ~exist('dbg','var')
    dbg = 0;
end

% Count control points in the root
nPoints = checkQuad(Boundary',NURBS.controlPoints');

% Split the root if there is more than one point in the root
if nPoints > 1
    % Set tree
    Quadtree = tree('root');
    data = Quadtree.Node{1,1};
    data = {data;...
           [];...
           NURBS};
      
    Quadtree = Quadtree.set(1, data);
    
    l = 1; % l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initialize at 1
    k = 0; % k: level of decomposition, equal to 0 at the root
    pos_aux = []; Q_aux = [0]; %auxiliar arrays
    
    [Quadtree] = decompose(Quadtree,Boundary,NURBS,l,k,pos_aux,Q_aux,k_min);
end

[Quadtree] = Star_Shape(Quadtree,NURBS);

[Quadtree] = QuadtreeBalance(Quadtree,NURBS);

% Balanced quadtree grid
% check if the quads intersected by the NURBS curve are star-shaped
% All quads intersected by the NURBS curve should be star-shaped !!
% Plot the kernels of quads intersected by the NURBS curve

if dbg == 1
    [Quadtree] = Star_Shape(Quadtree,NURBS,[],dbg);
end

end
