function [Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary)



% Count of control points in the root
nPoints = checkQuad( [min(Boundary(1,:)) min(Boundary(2,:))], [max(Boundary(1,:)) max(Boundary(2,:))],NURBS.controlPoints);

% Split the root if there is more than one point in the root
if nPoints > 1
    %Setting tree
    Quadtree = tree('root');
    data = Quadtree.Node{1,1};
    data = {data;...
          [];...
          NURBS.degree;...
          NURBS.knots;...
          NURBS.controlPoints;...
          NURBS.weights};
      
    Quadtree = Quadtree.set(1, data);
    
    l = 1; % l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initialize at 1
    k = 0; % k: level of decomposition, equal to 0 at the root
    pos_aux = []; Q_aux = [0]; %auxiliar arrays
    
    [Quadtree] = decompose(Quadtree,Boundary,NURBS,l,k,pos_aux,Q_aux,Boundary,k_min);
end

[Quadtree] = Star_Shape(Quadtree,NURBS,Boundary);

[Quadtree] = QuadtreeBalance(Quadtree,NURBS,Boundary);

end
