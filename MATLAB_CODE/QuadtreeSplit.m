function [Quadtree] = QuadtreeSplit(Quadtree,NURBS,seedingPoints, QuadLeaf_splitt)

%% cell decomposition
% number of seeding points
nSeedingPoints = size(seedingPoints,1);
splitted = [];
fatherWasSplitted =  0;

for i = 1:nSeedingPoints
    
    idx = 1;
    fatherWasSplitted =  0;
    while true
        
        idxChildren = Quadtree.getchildren(idx);
        
        if isempty(idxChildren)
            break
        end
        
        for j = 1:4
            
            Quad = Quadtree.Node{idxChildren(j),1}{7,1}(1:2,1:4);
                        
            if isPointInPolygon(Quad', seedingPoints(i,:))  ~= -1
                 
                idx = idxChildren(j);
                break
            end
        end
        
    end
  
    idxFather = Quadtree.Parent(idx);
    
    if  any (splitted == idxFather)
        fatherWasSplitted =  1;
    end
    
    if fatherWasSplitted ~= 1
        [Quadtree] = Decompose_helper(Quadtree,NURBS,idx);
        splitted = [splitted idx];
    end
    
end

%% local h-refinement of NURBS segment

numleaves = size(QuadLeaf_splitt,1);

for i = 1:numleaves
    idxLeaf = QuadLeaf_splitt(i,1);
    
    % get data stored in current leaf
    
    % 1. Quad name
    % 2. Quad location
    % 3. intersections points (physical space)
    % 4.
    % 5. intersection in parametric space of the curve
    % 6. NURBS definition
    %    NURBS degree
    %    NURBS control points
    %    NURBS Knot vector
    %    NURBS weights
    % 7. Quad definition
    % 8. Pointer to the children
    
    data = Quadtree.Node{idxLeaf,1};

%    % quad vertices
%     quadVertices = data{7};
%     % number of quad vertices
%     numQuadVertices = size(quadVertices,2);
%     % intersecion points
%     intersectionPoints = data{3};
%     % number of intersection points
%     numIntersectionPoints = size(intersectionPoints,2);
    
    if isempty(data{6})
        continue
    end

    % NURBS segment
    Segment = data{6,1};
    newSegment = hRefinement1d(Segment,1);
    data{6,1} = newSegment;
    Quadtree = Quadtree.set(idxLeaf,data);
end



end