l=Quadtree.findleaves();
function plot_leaf(Quadtree)
% plot_leaf: plots each quadtree's leaf

figure(3)
axis square
hold on
for i=1:length(l)
    
    intersections=Quadtree.Node{l(i),1}{5,1};
    % if empty plot empty quadrant
    if isempty(intersections) || length(intersections) == 1
         plot(Quadtree.Node{l(i),1}{10,1}(1,:),Quadtree.Node{l(i),1}{10,1}(2,:),'k.-');
    else
        % if each leaf may be plotted separately, uncomment
        %figure(i+1)
        
        % Obtain NURBS definition of the curve contained in each leaf
        NURBS.degree = Quadtree.Node{leaves(i),1}{6,1};
        NURBS.knots  = Quadtree.Node{leaves(i),1}{8,1};
        NURBS.controlPoints = Quadtree.Node{leaves(i),1}{7,1};
        NURBS.weights = Quadtree.Node{leaves(i),1}{9,1};
       
        NURBS_pts = CalculateNURBS(NURBS);
        plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2.5);
        
        % setting axis from function NURBS_parameter
%         axis(ax);
        
        
        % Plotting  leaf's quad
        plot(Quadtree.Node{l(i),1}{10,1}(1,:),Quadtree.Node{l(i),1}{10,1}(2,:),'k.-');
        
        % Plotting control points and control polygon for each leaf
        plot(Quadtree.Node{l(i),1}{7,1}(1,:),Quadtree.Node{l(i),1}{7,1}(2,:),'b-.','LineWidth',1)
        plot(Quadtree.Node{l(i),1}{7,1}(1,:),Quadtree.Node{l(i),1}{7,1}(2,:),'bo','LineWidth',1.5)
        
    end
end
end


