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



    end
    end
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


