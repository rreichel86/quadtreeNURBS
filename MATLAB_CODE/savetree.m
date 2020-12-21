function [Quadtree] = savetree(Q_aux, Quadtree,k, degree, Px, Py, ...
    newKnots, controlPoints, knots_new, weights,QS)
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


% Initializating variables 
controlPoints2=[];
knots_new2=[];
weights2=[];

% First level of decomposition
if k == 1
    for j=1:length(newKnotVals)/2
        % Loop over the inserted knots and extact NURBS segment 
        % contained in current quad
        k0 = find( abs(newKnots - newKnotVals(2*(j-1)+1)) < 1e-10, 1, 'First');
        k1 = find( abs(newKnots - newKnotVals(2*(j-1)+2)) < 1e-10, 1, 'Last');
        
        if k0 == 1 && k1 ~= length(newKnots)
            controlPoints2 = controlPoints(1:2,k0:(k1-degree));
            knots_new2 = [newKnots(k0:k1) newKnots(k1)];
            weights2 = weights(k0:(k1-degree));
        elseif k0 ~= 1 && k1 == length(newKnots)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree-1));
            knots_new2 = [newKnots(k0) newKnots(k0:k1)];
            weights2 = weights((k0-1):(k1-degree-1));
        elseif k0~=1 && k1 ~= length(newKnots)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree));
            knots_new2 = [newKnots(k0) newKnots(k0:k1) newKnots(k1)];
            weights2 = weights((k0-1):(k1-degree));
        else    
            controlPoints2 = controlPoints;
            knots_new2 = newKnots;
            weights2 = weights;
        end
    end
    
    % Store information in tree's new node 
    
    % Quad name 
    % Quad location
    % horizontal intersections (physical space)
    % vertical intersections (physical space)
    % intersection in parametric space of the curve
    % NURBS degree
    % NURBS control points 
    % NURBS Knot vector 
    % NURBS weights
    % Quad definition 
    % Pointer to the children
    data = {['Quad' num2str(Q_aux)];...
             Q_aux;...
             Px;...
             Py;...
             newKnotVals;...
             degree;...
             controlPoints2;...
             knots_new2;...
             weights2;...
             QS;...
             []};
         
    % attach tree's new node to the root     
    [Quadtree, nodeID] = Quadtree.addnode(1, data);
    % add nodeID to root's pointer list 
    data = Quadtree.Node{1,1};
    data{2} = [data{2} nodeID];
    Quadtree = Quadtree.set(1, data);
else
    % other levels of decomposition
    for j=1:length(newKnotVals)/2
        % Loop over the inserted knots and extact NURBS segment 
        % contained in current quad
        k0 = find( abs(newKnots - newKnotVals(2*(j-1)+1)) < 1e-10, 1, 'First');
        k1 = find( abs(newKnots - newKnotVals(2*(j-1)+2)) < 1e-10, 1, 'Last');
        
        
        if k0 == 1 && k1 ~= length(newKnots)
            controlPoints2 = controlPoints(1:2,k0:(k1-degree));
            knots_new2 = [newKnots(k0:k1) newKnots(k1)];
            weights2 = weights(k0:(k1-degree));
        elseif  k0 ~= 1 && k1 == length(newKnots)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree-1));
            knots_new2 = [newKnots(k0) newKnots(k0:k1)];
            weights2 = weights((k0-1):(k1-degree-1));
        elseif k0 ~= 1 && k1 ~= length(newKnots)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree));
            knots_new2 = [newKnots(k0) newKnots(k0:k1) newKnots(k1)];
            weights2 = weights((k0-1):(k1-degree));
        else    
            controlPoints2 = controlPoints;
            knots_new2 = newKnots;
            weights2 = weights;    
        end
    end
    
    % Store information in tree's new node 
    
    % Quad name 
    % Quad location
    % horizontal intersections (physical space)
    % vertical intersections (physical space)
    % intersection in parametric space of the curve
    % NURBS degree
    % NURBS control points 
    % NURBS Knot vector 
    % NURBS weights
    % Quad definition 
    % Pointer to the children
    data = {['Quad' num2str(Q_aux)];...
            Q_aux;...
            Px;...
            Py;...
            newKnotVals;...
            degree;...
            controlPoints2;...
            knots_new2;...
            weights2;...
            QS;...
            []};
           
    refFQ = Q_aux(1:end-2);
    Parent = findQ(Quadtree,refFQ);
    % attach tree's new node to its father
    [Quadtree, nodeID] = Quadtree.addnode(Parent, data);
    % add nodeID to father's pointer list 
    data = Quadtree.Node{Parent,1};
    data{11} = [data{11} nodeID];
    Quadtree = Quadtree.set(Parent, data);
    
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
