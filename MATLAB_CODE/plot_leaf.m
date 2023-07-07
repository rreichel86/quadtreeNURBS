function plot_leaf(Quadtree,plotOptions)
% plot_leaf: plots each quadtree's leaf

if ~exist('plotOptions','var')
    % 1 - plot NURBS curve 
    % 2 - plot Quad 
    % 3 - plot control polygon and control points
    % 4 - plot intersection points
    plotOptions = [0 1 0 0];
end


% Get Quadtree leaves
leaves = Quadtree.findleaves();
numLeaves = length(leaves);

figure(2)
% xticks([])
% yticks([])
daspect([1 1 1])
box on
set(gca,'TickLabelInterpreter','latex','FontSize',20,'FontName','Times');
axis square
hold on

for i = 1: numLeaves
    
    idxLeaf = leaves(i);
    
    % get data stored in current leaf
    
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
    
    data = Quadtree.Node{idxLeaf,1};
    
    intersections = data{5};
    Quad = data{7};
    
    % check if leaf has intersections
    if isempty(intersections)
        patch(Quad(1,:),Quad(2,:),'w','FaceAlpha',0,'LineStyle','-','LineWidth',1);
    elseif length(intersections) == 1  
        patch(Quad(1,:),Quad(2,:),'r','FaceAlpha',0,'LineStyle','-','LineWidth',1);
    else

        % Obtain definition of the NURBS curve contained in each leaf
        NURBS = data{6};

        ncp = length(NURBS.knots) - NURBS.degree - 1;
        NURBS_pts = CalculateNURBS(NURBS);

        % Plot NURBS curve
        if plotOptions(1) == 1
            plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2);
        end
        % Plot Quad
        if plotOptions(2) == 1
            patch(Quad(1,:),Quad(2,:),'w','FaceAlpha',0,'LineStyle','-','LineWidth',1);
        end
        % Plot control polygon
        if plotOptions(3) == 1
            plot(NURBS.controlPoints(1,:),NURBS.controlPoints(2,:),'b-.','LineWidth',1)
            % Plot control points
%             plot(NURBS.controlPoints(1,2:ncp-1),NURBS.controlPoints(2,2:ncp-1),'o','Color','red','MarkerFaceColor','r','MarkerSize',6)
            plot(NURBS.controlPoints(1,1:ncp),NURBS.controlPoints(2,1:ncp),'o','Color','red','MarkerFaceColor','r','MarkerSize',6)
        end
        % Plot control points/ plot intersection points
        if plotOptions(4) == 1
            plot(NURBS.controlPoints(1,1),NURBS.controlPoints(2,1),'bo','LineWidth',1.5)
        end

    end
end
% hold off

end


