function [Quadtree] = decompose(Quadtree, Quad,  controlPoints, knots,...
    weights, degree, l, k,pos_aux, Q_aux, Boundary)
% decompose function divides a given quad into 4 children. It calls the
% splitting function for obtaining the possible NURBS segment contained
% in each children.
%
% INPUT:
% Quadtree data
% Quad: geometrical definition of the quad that is about to be splitted
% Definition of the NURBS
% l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initializate at 1
% k: level of decomposition
%
% OUTPUT:
% Quadtree data after performing decomposition and eventual curve splitting


% Each time decompose is called the var k is increased, since we are splitting
% the quad
k=k+1;


% Split the parent quad into 4 child quads
[child1 , child2 , child3 , child4] = splitQuad(min(Quad(1,:)),...
    max(Quad(1,:)), min(Quad(2,:)), max(Quad(2,:)) );
child = [child1 , child2 , child3 , child4];


for i = 1:4
    % Loop over the 4 children. First we split the eventual NURBS segment
    % contained in the given children. Then we check if there is more than
    % one seed point or 2 segments of the curve in the quad. If one of
    % those two conditions is true, we  decompose the given quad
    
    % Selecting given quad
    Current = [child(1,2*i-1),child(1,2*i);child(2,2*i-1),child(2,2*i)];
    Quad = [Current(1,1),Current(1,2),Current(1,2),Current(1,1),Current(1,1);...
        Current(2,1),Current(2,1),Current(2,2),Current(2,2),Current(2,1)];
    plot([Quad(1,:), Quad(1,1)],[Quad(2,:), Quad(2,1)],'k')
    
    %Call the splitting function. Get the tree data structure after
    %performing the splitting at the given quad and an auxiliary array
    [Q_aux, Quadtree,numInterPoints] = splitting(Quadtree, Q_aux, Quad, controlPoints,...
        knots, weights, degree,l, k, i, pos_aux);
 
    %Obtains number of seed points in given quad
    nPointsCurrent = checkQuad(Current(:,1)', Current(:,2)', controlPoints);
    
    %If there are more than two segments or seed points, decompose again
    if nPointsCurrent > 1 || numInterPoints > 2
        
        [Quadtree] = decompose(Quadtree,Quad, controlPoints, knots,...
            weights, degree, l, k,pos_aux, Q_aux,Boundary);
    end
    
    if k == 1
        l=l+1;
    end
end
end
