function plot_leaf(Quadtree,ax)
% Function plot_leaf plots each quadtree's leaf content
l=Quadtree.findleaves();
figure(3)
hold on
for i=1:length(l)
    
    intersections=Quadtree.Node{l(i),1}{5,1};
    % if empty plot empty quadrant
    if isempty(intersections) || length(intersections) == 1
%         plot(Quadtree.Node{l(i),1}{10,1}(1,:),Quadtree.Node{l(i),1}{10,1}(2,:),'k');
    else
        % if each leaf may be plotted separately, uncomment
        %figure(i+1)
        
        
        % setting axis from function NURBS_parameter
%         axis(ax);
        
        % Obtaining NURBS definition of the curve contained in each leaf quad
        leave=CalculateNURBS( Quadtree.Node{l(i),1}{6,1},Quadtree.Node{l(i),1}{8,1}, Quadtree.Node{l(i),1}{7,1},  Quadtree.Node{l(i),1}{9,1});
        plot(leave(:,1),leave(:,2),'r','LineWidth',3);
        
        % Plotting  leaf's quad
        plot(Quadtree.Node{l(i),1}{10,1}(1,:),Quadtree.Node{l(i),1}{10,1}(2,:),'k');
        
        % Plotting control points and control polygon for each leaf
        plot(Quadtree.Node{l(i),1}{7,1}(1,:),Quadtree.Node{l(i),1}{7,1}(2,:),'--','LineWidth',1)
        plot(Quadtree.Node{l(i),1}{7,1}(1,:),Quadtree.Node{l(i),1}{7,1}(2,:),'bo','LineWidth',1)
        
    end
end
end


