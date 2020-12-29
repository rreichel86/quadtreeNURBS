function [Quadtree] = savetree(Q_aux, Quadtree,k, Px, Py, knotVals,...
    NURBS, QS)
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

% NURBS segment
NURBS_segment = struct;

if ~isempty(NURBS)
    
    degree = NURBS.degree;
    knots = NURBS.knots;
    controlPoints = NURBS.controlPoints;
    weights = NURBS.weights;
    
    % Loop over the inserted knots and extact NURBS segment
    % contained in current quad
    for j = 1:length(knotVals)/2
        
        k0 = find( abs(knots - knotVals(2*(j-1)+1)) < 1e-10, 1, 'First');
        k1 = find( abs(knots - knotVals(2*(j-1)+2)) < 1e-10, 1, 'Last');
        
        if k0 == 1 && k1 ~= length(knots)
            NURBS_segment.degree = degree;
            NURBS_segment.controlPoints = controlPoints(1:2,k0:(k1-degree));
            NURBS_segment.knots = [knots(k0:k1) knots(k1)];
            NURBS_segment.weights = weights(k0:(k1-degree));
        elseif k0 ~= 1 && k1 == length(knots)
            NURBS_segment.degree = degree;
            NURBS_segment.controlPoints = controlPoints(1:2,(k0-1):(k1-degree-1));
            NURBS_segment.knots = [knots(k0) knots(k0:k1)];
            NURBS_segment.weights = weights((k0-1):(k1-degree-1));
        elseif k0 ~= 1 && k1 ~= length(knots)
            NURBS_segment.degree = degree;
            NURBS_segment.controlPoints = controlPoints(1:2,(k0-1):(k1-degree));
            NURBS_segment.knots = [knots(k0) knots(k0:k1) knots(k1)];
            NURBS_segment.weights = weights((k0-1):(k1-degree));
        else
            NURBS_segment.degree = degree;
            NURBS_segment.controlPoints = controlPoints;
            NURBS_segment.knots = knots;
            NURBS_segment.weights = weights;
        end
    end
end
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
    Quadtree = Quadtree.set(Parent, data);   
end

% Deleting information of the father, information contained in leafs
%     if length(Quadtree.Node{Parent,1}{11}) == 4
%         data = Quadtree.Node{Parent,1};
%         data{6} = [];
%         data{7} = [];
%         data{8} = [];
%         data{9} = [];
%         Quadtree = Quadtree.set(Parent, data);
%     end

end
