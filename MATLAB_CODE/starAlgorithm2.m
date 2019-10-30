function [Quadtree]=starAlgorithm2(Quadtree,f)

% starAlgorithm2 finds kernels of a quadtree leaf
% using the coresponding polygonised NURBS curve in the quad
% for splitting into 2 subdomains

% Input: 
% Quadtree 
% f: current leaf to be searched for kernel

% Output: 
% Quadtree: updated Quadtree, containing the coordinates of vertices of the
% kernels at rows 12 and 13 for each leaf node in Quadtree

% Calculate NURBS curve in current quad
curve=CalculateNURBS( Quadtree.Node{f,1}{6,1}, Quadtree.Node{f,1}{8,1}, Quadtree.Node{f,1}{7,1}, Quadtree.Node{f,1}{9,1}); % obtain NURBS from current Quad

%selects n points from the curve, taking into account the degree of the
%NURBS curve. Ex.: for curve of degree 3, n = 3 + 10 = 13
n_points=floor(linspace(1,length(curve(1,:)),Quadtree.Node{f,1}{6,1}+10));

% Constructing the polygonised curve with n coordinates of the curve
% In reallity it is a polygon, but with close approximation to the real
% NURBS curve
for j=1:length(n_points)
    points(:,j)=curve(1:2,n_points(j));
end

points=round(points,5);

% Obtain intersections of control polygon with quad from data tree
Int = [Quadtree.Node{f,1}{3,1}, Quadtree.Node{f,1}{4,1}];
Int = round(Int,5);

% Store coordinates of current quad in variable
Quad_current = Quadtree.Node{f,1}{10,1}(:,1:length(Quadtree.Node{f,1}{10,1})-1);
Quad_current = round(Quad_current,5);

% Call the function starSplitPolygons, which will obtain the 2 subdomains
% after their separation by the polygonised curve
[polygon1, polygon2] = starSplitPolygons(Quad_current,points,Int);

% Plot Leaves and polygonised curve
% figure(7);
% hold on
% plot(points(1,:),points(2,:),'black')
% Quad_current_closed = [Quad_current, Quad_current(:,1)];
% plot(Quad_current_closed(1,:),Quad_current_closed(2,:),'black','LineWidth', 1.5)
% title('Leaves and polygonised curve')


% PLot polygons 1 and 2 in different colors
% polygon1_closed = [polygon1, polygon1(:,1)];
% polygon2_closed = [polygon2, polygon2(:,1)];
% figure(8);
% hold on
% plot(polygon1(1,:),polygon1(2,:),'red')
% plot(polygon2(1,:),polygon2(2,:),'green')
% title('Polygon 1 and 2')

% Runnung the kernel search for polygon 1 & 2
K1 = starShapedCheck(polygon1);
K2 = starShapedCheck(polygon2);

% Plot the kernels 1 and 2 in different colors
% figure(9)
% hold on;
% if ~isempty(K1)
%     K1_closed = [K1, K1(:,1)];
%     plot(K1_closed(1,:),K1_closed(2,:),'black')
% end
% if ~isempty(K2)
%     K2_closed = [K2, K2(:,1)];
%     plot(K2_closed(1,:),K2_closed(2,:),'black')
% end
% title('K1 and K2')

% Store K1 and K2 in Quadtree
data=Quadtree.Node{f,1};
data{12}=K1;
data{13}=K2;
Quadtree = Quadtree.set(f, data);

end





