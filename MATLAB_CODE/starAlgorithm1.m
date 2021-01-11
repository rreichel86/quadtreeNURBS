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

% get data stored in current leaf

% 1. Quad name
% 2. Quad location
% 3. intersection points (physical space)
% 4.
% 5. intersection in parametric space of the curve
% 6. NURBS definition
%    NURBS degree
%    NURBS control points
%    NURBS Knot vector
%    NURBS weights
% 7. Quad definition
% 8. Pointer to the children

data = Quadtree.Node{idx,1};

% quad vertices 
quadVertices = data{7};
% number of quad vertices
numQuadVertices = size(quadVertices,2);
% control points
NURBS_segment = data{6};
controlPoints = NURBS_segment.controlPoints; 
% intersecion points
intersectionPoints = data{3}; 
% number of intersection points
numIntersectionPoints = size(intersectionPoints,2);

% number of polygon's vertices
numPolyVertices = numQuadVertices + numIntersectionPoints;
% polygon vertices
poly = zeros(2, numPolyVertices); 
% intersection points location 
locIntersectionPoints = zeros(1,2);

nPv = 0;
for iQv = 1:numQuadVertices 
    
    A = quadVertices(:,iQv);
    
    nPv = nPv + 1;
    poly(:,nPv) = A;
    
    if iQv ~= numQuadVertices
        B = quadVertices(:,iQv+1);
    else
        B = quadVertices(:,1);
    end
    
    for iIp = 1:numIntersectionPoints
        
        if locIntersectionPoints(iIp) ~= 0
           continue 
        end    
        
        P = intersectionPoints(1:2,iIp);
        
        [pointInSegment, locPoint] = isPointInLineSegment(A,B,P,3);
        
        if pointInSegment
            
            if locPoint > 0
                nPv = nPv + 1;
                poly(:,nPv) = P;
                locIntersectionPoints(iIp) = nPv;
            else 
                locIntersectionPoints(iIp) = nPv;
            end   
            
        end  
    end
end

% update number of polygon's vertices
if (nPv ~= numPolyVertices)
    numPolyVertices = nPv;
end     

% insert control points and split quad 
I1 = poly(:,locIntersectionPoints(1));

if norm( I1 - controlPoints(:,1) ) < tol 
    polygon_1 = [poly(:,1:locIntersectionPoints(1)),...
                 controlPoints(:,2:end-1),...
                 poly(:,locIntersectionPoints(2):numPolyVertices)];
             
    polygon_2 = [poly(:,locIntersectionPoints(1):locIntersectionPoints(2)),...
                 controlPoints(:,end-1:-1:1)];
else
    polygon_1 = [poly(:,1:locIntersectionPoints(1)),...
                 controlPoints(:,end-1:-1:2),...
                 poly(:,locIntersectionPoints(2):numPolyVertices)];
             
    polygon_2 = [poly(:,locIntersectionPoints(1):locIntersectionPoints(2)),...
                 controlPoints(:,2:end)];
end

% compute polygons kernel
K1 = computePolygonKernel(polygon_1');
K2 = computePolygonKernel(polygon_2');
star_shaped = ~isempty(K1) && ~isempty(K2);
end


