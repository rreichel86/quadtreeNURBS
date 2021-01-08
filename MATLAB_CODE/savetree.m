function [Quadtree] = savetree(Q_aux, Quadtree,k, Px, Py, knotVals,...
    NURBS_segment, QS)
% Savetree function prepares the information that will be store at the tree
% data structure. After splitting the function the description of the NURBS
% contained in the quad is stored in a node of the tree data structure,
% with it respectives pointers and auxiliar data. The code from Jean-Yves
% Tinevez is used for generating a tree data as a matlab class.
%
% Input:
% Quadtree data
% Quad: geometrical definition of the considered quad
% Definition of the NURBS
% Auxiliar variables previously defined
% Output:
% Quadtree after storing new information


% Store information in tree's new node

% 1. Quad name
% 2. Quad location
% 3. horizontal intersections (physical space)
% 4. vertical intersections (physical space)
% 5. intersection in parametric space of the curve
% 6. NURBS definition
%    NURBS degree
%    NURBS control points
%    NURBS Knot vector
%    NURBS weights
% 7. Quad definition
% 8. Pointer to the children
data = {['Quad ' num2str(Q_aux)];...
    Q_aux;...
    Px;...
    Py;...
    knotVals;...
    NURBS_segment;...
    QS;...
    []};


if k == 1 % first level 
    
    % attach tree's new node to the root
    [Quadtree, nodeID] = Quadtree.addnode(1, data);
    % add nodeID to root's pointer list
    data = Quadtree.Node{1,1};
    data{2} = [data{2} nodeID];
    Quadtree = Quadtree.set(1, data);
else % other levels 
    
    refFQ = Q_aux(1:end-2);
    Parent = findQ(Quadtree,refFQ);
    % attach tree's new node to its father
    [Quadtree, nodeID] = Quadtree.addnode(Parent, data);
    % add nodeID to father's pointer list
    data = Quadtree.Node{Parent,1};
    data{8} = [data{8} nodeID];
     
    % delete stored information in the father
    
    % 3. horizontal intersections (physical space)
    % 4. vertical intersections (physical space)
    % 5. intersection in parametric space of the curve
    % 6. NURBS definition
    if length(data{8}) == 4
        data{3} = [];
        data{4} = [];
        data{5} = [];
        data{6} = [];
    end    
    
    Quadtree = Quadtree.set(Parent, data); 
    
end


end
