function [Quadtree,polygon1, polygon2]=starAlgorithm1(Quadtree,f)

% starAlgorithm1 finds kernels of a quadtree leaf
% using the coresponding control polygon in the quad
% for splitting into 2 subdomains

% Input: 
% Quadtree 
% f: current leaf to be searched for kernel

% Output: 
% Quadtree: updated Quadtree, containing the coordinates of vertices of the
% kernels at rows 12 and 13 for each leaf node in Quadtree

% first obtain the control polygon and store in variable "points"
points=round(Quadtree.Node{f,1}{7,1},5); % Store CPs in variable

% Obtain intersections of control polygon with quad from data tree
Int = [Quadtree.Node{f,1}{3,1}, Quadtree.Node{f,1}{4,1}];
Int = round(Int,5);

% Store coordinates of current quad in variable
Quad_current = Quadtree.Node{f,1}{10,1}(:,1:length(Quadtree.Node{f,1}{10,1})-1);
Quad_current = round(Quad_current,5);

% Call the function starSplitPolygons, which will obtain the 2 subdomains
% after their separation by the control polygon
[polygon1, polygon2] = starSplitPolygons(Quad_current,points,Int);


% Plot Leaves and control polygon
%  figure(7);
%  hold on
%  plot(points(1,:),points(2,:),'black')
%  Quad_current_closed = [Quad_current, Quad_current(:,1)];
%  plot(Quad_current_closed(1,:),Quad_current_closed(2,:),'black', 'LineWidth', 1.5)
%  title('Leaves and control polygon')

 %PLot polygons 1 and 2 in different colors
 polygon1_closed = [polygon1, polygon1(:,1)];
 polygon2_closed = [polygon2, polygon2(:,1)];
 figure(8);
 hold on
 plot(polygon1_closed(1,:),polygon1_closed(2,:),'red')
 plot(polygon2_closed(1,:),polygon2_closed(2,:),'green')
 title('Polygon 1 and 2')

% Runnung the kernel search for polygon 1 & 2
K1 = starShapedCheck(polygon1);
K2 = starShapedCheck(polygon2);

% Plot the kernels 1 and 2 in different colors
% figure(9)
% hold on;
% color_fill = [1.0000    0.3608    0.3608];
% if ~isempty(K1)
%     K1_closed = [K1, K1(:,1)];
%     fill(K1_closed(1,:),K1_closed(2,:),'red')
% end
% if ~isempty(K2)
%     K2_closed = [K2, K2(:,1)];
%     fill(K2_closed(1,:),K2_closed(2,:),'green')
% end
% title('K1 and K2')
%

% Store K1 and K2 in Quadtree
data=Quadtree.Node{f,1};
data{12}=K1;
data{13}=K2;
Quadtree = Quadtree.set(f, data);

% elements(Quadtree,polygon1, polygon2)

end




