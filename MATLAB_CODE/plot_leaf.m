function plot_leaf(Quadtree)
% plot_leaf: plots each quadtree's leaf

% Get Quadtree leaves
leaves = Quadtree.findleaves();
numLeaves = length(leaves);

figure(3)
axis square
hold on

for i = 1:numLeaves
    
    intersections=Quadtree.Node{leaves(i),1}{5,1};
    
    % check if leaf has intersections
    if isempty(intersections) || length(intersections) == 1
        % Plot Quad
        plot(Quadtree.Node{leaves(i),1}{10,1}(1,:),Quadtree.Node{leaves(i),1}{10,1}(2,:),'k.-');
    else
        % Obtain NURBS definition of the curve contained in each leaf
        NURBS.degree = Quadtree.Node{leaves(i),1}{6,1};
        NURBS.knots  = Quadtree.Node{leaves(i),1}{8,1};
        NURBS.controlPoints = Quadtree.Node{leaves(i),1}{7,1};
        NURBS.weights = Quadtree.Node{leaves(i),1}{9,1};
       
        NURBS_pts = CalculateNURBS(NURBS);
        plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2.5);
        
        % Plot Quad
        plot(Quadtree.Node{leaves(i),1}{10,1}(1,:),Quadtree.Node{leaves(i),1}{10,1}(2,:),'k.-');
        
        % Plot control points and control polygon for each leaf
        plot(NURBS.controlPoints(1,:),NURBS.controlPoints(2,:),'b-.','LineWidth',1)
        plot(NURBS.controlPoints(1,:),NURBS.controlPoints(2,:),'bo','LineWidth',1.5)
        
    end
end
hold off

end


