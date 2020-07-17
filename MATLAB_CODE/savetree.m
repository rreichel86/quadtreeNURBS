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

% Storing data for the first level of decomposition
if k==1
    for j=1:length(newKnots)/2
        % Loop over the inserted knots and selecting section of the NURBS
        % contained in the quad
        k0 = find( abs(knots_new - newKnots(2*(j-1)+1)) < 1e-10, 1, 'First');
        k1 = find( abs(knots_new - newKnots(2*(j-1)+2)) < 1e-10, 1, 'Last');
        
        if k0==1
            controlPoints2 = controlPoints(1:2,k0:(k1-degree));
            knots_new2= [knots_new(k0:k1) knots_new(k1)];
            weights2= weights(k0:(k1-degree));
        end
        if k1==length(knots_new)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree-1));
            knots_new2= [knots_new(k0) knots_new(k0:k1)];
            weights2= weights((k0-1):(k1-degree-1));
        end
        if k0~=1 && k1~=length(knots_new)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree));
            knots_new2= [knots_new(k0) knots_new(k0:k1) knots_new(k1)];
            weights2= weights((k0-1):(k1-degree));
        end
    end
    % Storing information at tree's new node and attaching to root
    data={['Quad' num2str(Q_aux)];Q_aux; Px; Py; newKnots; degree; ...
        controlPoints2; knots_new2; weights2;QS;[]};
    [Quadtree, node1] = Quadtree.addnode(1, data);
    % Adding pointer to root's new children
    data=Quadtree.Node{1,1};
    data{2}=[data{2} node1];
    Quadtree = Quadtree.set(1, data);
else
    % Storing data for any other level of decomposition
    for j=1:length(newKnots)/2
        % Loop over the inserted knots and selecting section of the NURBS
        % contained in the quad
        k0 = find( abs(knots_new - newKnots(2*(j-1)+1)) < 1e-10, 1, 'First');
        k1 = find( abs(knots_new - newKnots(2*(j-1)+2)) < 1e-10, 1, 'Last');
        
        
        if k0==1
            controlPoints2 = controlPoints(1:2,k0:(k1-degree));
            knots_new2= [knots_new(k0:k1) knots_new(k1)];
            weights2= weights(k0:(k1-degree));
        end
        if k1==length(knots_new)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree-1));
            knots_new2= [knots_new(k0) knots_new(k0:k1)];
            weights2= weights((k0-1):(k1-degree-1));
        end
        if k0~=1 && k1~=length(knots_new)
            controlPoints2 = controlPoints(1:2,(k0-1):(k1-degree));
            knots_new2= [knots_new(k0) knots_new(k0:k1) knots_new(k1)];
            weights2= weights((k0-1):(k1-degree));
        end
    end
    
    % Storing information at tree's new node and attaching to father
    data={['Quad' num2str(Q_aux)];Q_aux; Px ;Py; newKnots; degree;...
        controlPoints2; knots_new2; weights2;QS;[]};
    for i=1:length(Quadtree.Node);Pos{i}=Quadtree.Node{i,1}{2};end
    Q_aux=Q_aux(1:end-2);
    idx = cellfun('length',Pos)==length(Q_aux);
    for j=1:length(idx);tf=isequal(Pos{j},Q_aux);if tf;Parent=j;end;end
    [Quadtree, node2] = Quadtree.addnode(Parent, data);
    
    % Adding pointer to father's new child
    data=Quadtree.Node{Parent,1};
    data{11}=[data{11} node2];
    Quadtree = Quadtree.set(Parent, data);
    
    % Deleting information of the father, information contained in leafs
    if length(Quadtree.Node{Parent,1}{11})==4
    data=Quadtree.Node{Parent,1};
    data{6}=[];data{7}=[];data{8}=[];data{9}=[];
    Quadtree = Quadtree.set(Parent, data);      
    end
    
end
