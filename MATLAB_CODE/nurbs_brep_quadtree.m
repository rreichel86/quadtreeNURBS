function [Quadtree] = nurbs_brep_quadtree(k_min,degree,knots,controlPoints,weights,Boundary)



% Count of control points in the root
nPoints = checkQuad( Boundary(:,2)',  Boundary(:,4)',controlPoints);

% Split the root if there is more than one point in the root
if nPoints > 1
    %Setting tree
    Quadtree = tree('root');
    data=Quadtree.Node{1,1};
    data={data;[];degree;knots;controlPoints;weights};
    Quadtree = Quadtree.set(1, data);
    
    l=1; % l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initialize at 1
    k=0; % k: level of decomposition, equal to 0 at the root
    % k_min = 2;
    pos_aux=[]; Q_aux=[0]; %auxiliar arrays
    
    [Quadtree] = decompose(Quadtree,Boundary,controlPoints, knots,...
        weights, degree, l,k,pos_aux, Q_aux,Boundary,k_min);
end

[Quadtree] = Star_Shape(Quadtree,controlPoints, knots, weights,...
    degree, Boundary);

[Quadtree] = QuadtreeBalance(Quadtree, controlPoints, knots, ...
    weights, degree,Boundary);
% [Quadtree] = Balance_Quadtree(Quadtree, controlPoints, knots, ...
%     weights, degree,Boundary);


end
