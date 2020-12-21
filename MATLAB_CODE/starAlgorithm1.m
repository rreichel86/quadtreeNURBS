function [star_shaped, Quadtree] = starAlgorithm1(Quadtree,idx)
% starAlgorithm1: compute kernels of a Quadtree leaf
% Only for Quadtree leaves that are splitted by the NURBS curve
%
% INPUT: 
% Quadtree ------------------- Quadtree data structure
% idx ------------------------ current leaf index
%
% OUTPUT: 
% Quadtree ------------------- updated Quadtree data structure
% star_shaped ---------------- 1: if and only if both Quad subregions are 
%                                 star-shaped
%                              0: otherwise
%
% -------------------------------------------------------------------------

tol = 1e-10;
% Quad vertices
vertices = Quadtree.Node{idx,1}{7,1}(:,1:end)';
% control points
controlPoints = Quadtree.Node{idx,1}{6,1}.controlPoints';
% vertical intersections [left; right] 
vertIntrsc = Quadtree.Node{idx,1}{4,1}';
nVertIntrsc = size(vertIntrsc,1);
% horizontal intersections [bottom; top]
horzIntrsc = Quadtree.Node{idx,1}{3,1}';
nHorzIntrsc = size(horzIntrsc,1);
% intersections
intrsc = [horzIntrsc;vertIntrsc];

% Quad left corner
Qxmin = min(vertices(:,1));
Qymin = min(vertices(:,2));

% intersection points indices
iintrsc = zeros(2,1);
poly = zeros(9,2) - 99;
poly(1:2:8,:) = vertices;
poly(9,:) = vertices(1,:);

if nVertIntrsc == 2
    poly(4,:) = intrsc(2,:);
    poly(8,:) = intrsc(1,:); 
    iintrsc = [4,8];
elseif nHorzIntrsc == 2
    poly(2,:) = intrsc(1,:); 
    poly(6,:) = intrsc(2,:);
    iintrsc = [2,6];
else
    % horizontal intersection
    if abs(intrsc(1,2) - Qymin) < tol % bottom 
        poly(2,:) = intrsc(1,:);
        iintrsc(1) = 2;
    else % top
        poly(6,:) = intrsc(1,:);
        iintrsc(1) = 6;
    end
    % vertical intersection
    if abs(intrsc(2,1) - Qxmin) < tol % left 
        poly(8,:) = intrsc(2,:);
        iintrsc(2) = 8;
    else % right
        poly(4,:) = intrsc(2,:);
        iintrsc(2) = 4;
    end
end

% swap intersection points indices
if iintrsc(1) > iintrsc(2) 
    tmp = iintrsc(1);
    iintrsc(1) = iintrsc(2);
    iintrsc(2) = tmp;
end

% insert control point
if norm( poly(iintrsc(1),:) - controlPoints(1,:) ) < tol 
    polygon_1 = [poly(1:iintrsc(1),:);...
                 controlPoints(2:end-1,:);...
                 poly(iintrsc(2):end,:)];
             
    polygon_2 = [poly(iintrsc(1):iintrsc(2),:);...
                 controlPoints(end-1:-1:1,:);...
                ];
else
    polygon_1 = [poly(1:iintrsc(1),:);...
                 controlPoints(end-1:-1:2,:);...
                 poly(iintrsc(2):end,:)];
             
    polygon_2 = [poly(iintrsc(1):iintrsc(2),:);...
                 controlPoints(2:end,:);...
                ];
end

polygon_1 = polygon_1(polygon_1(:,1) ~= -99,:);
polygon_2 = polygon_2(polygon_2(:,1) ~= -99,:);

% Compute polygons kernel
K1 = computePolygonKernel(polygon_1);
K2 = computePolygonKernel(polygon_2);
star_shaped = ~isempty(K1) && ~isempty(K2);
end


