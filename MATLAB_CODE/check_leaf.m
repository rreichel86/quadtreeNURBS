function [Quadtree] = check_leaf(Quadtree,dbg)
% check_leaf


if ~exist('dbg','var')
    dbg = 0;
end

if dbg == 1
    figure(2)
    axis square
    hold on
end 

tol = 1e-10;
% get Quadtree leaves
leaves = Quadtree.findleaves();
numLeaves = length(leaves);

updatedLeaves = [];

for i = 1:numLeaves
    
    idxLeaf = leaves(i);
    
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
    
    data = Quadtree.Node{idxLeaf,1};
    
    
    refLeaf = data{2};
    levelLeaf = length(refLeaf);
    
    % intersections
    intersections = data{5};
    
    % check if leaf has intersections
    if isempty(intersections) || length(intersections) == 1
        continue
    else
        % NURBS segment
        NURBS = data{6};
        % quad vertices
        quadVertices = data{7};
        % number of quad vertices
        numQuadVertices = size(quadVertices,2);
        % intersecion points
        intersectionPoints = data{3};
        % number of intersection points
        numIntersectionPoints = size(intersectionPoints,2);
        % intersection points location
        intersectedEdges = data{4};
        alphas = zeros(1,2);
        
        if isequal( intersectedEdges, [1 4] )
            vertexNum = 1;
        elseif isequal( intersectedEdges, [1 2] )
            vertexNum = 2;
        elseif isequal( intersectedEdges, [2 3] )
            vertexNum = 3;
        elseif isequal( intersectedEdges, [3 4] )
            vertexNum = 4;
        else
            continue
        end
        
        for iIp = 1:numIntersectionPoints
            
            iQv = intersectedEdges(iIp);
            
            Vi = iQv;
            if iQv ~= numQuadVertices
                Vj = iQv+1;
            else
                Vj = 1;
            end
            
            if Vi == vertexNum
                A = quadVertices(:,Vi);
                B = quadVertices(:,Vj);
            else
                A = quadVertices(:,Vj);
                B = quadVertices(:,Vi);
            end
            
            P = intersectionPoints(1:2,iIp);
            
            [pointInSegment, locPoint] = isPointInLineSegment(A,B,P,4);
            
            if pointInSegment
                alphas(iIp) = locPoint;
            end
            
        end
        
        if alphas(1) < 1/5 && alphas(2) < 1/5
            
            % project quad vertex onto NURBS segment
            oldVertexCoor = quadVertices(:,vertexNum);
            u0 = intersectionPoints(3,1);
            
            [newVertexCoor,u] = pointProjectionNURBS(oldVertexCoor,u0,NURBS);
            
            % replace quad vertex by its projection
            quadVertices(:,vertexNum) = newVertexCoor;
            
            % plot newVertex and quad (current leaf)
            if dbg == 1
                plot(quadVertices(1,vertexNum),quadVertices(2,vertexNum), 'o','Color','k','MarkerFaceColor','k','MarkerSize',8);
                patch(quadVertices(1,:),quadVertices(2,:),'r','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
            end
            
            % update intersection points 
            data{3} = [newVertexCoor; u];
            data{4} = vertexNum;
            data{5} = u;
            data{6} = [];
            % update quad vertices 
            data{7} = quadVertices;
            
            % update data stored in current leaf
            Quadtree = Quadtree.set(idxLeaf, data);
            
            LEAF.idx = idxLeaf;
            LEAF.vertexNum = vertexNum;
            LEAF.vertexCoor = [newVertexCoor; u];
            LEAF.oldVertexCoor = oldVertexCoor;
            LEAF.oldIntersectionPoints = intersectionPoints;
            
            updatedLeaves = [updatedLeaves LEAF];
            
            % seach Neighbours that share the same vertex
            % Vertex 1
            %         10. West Neighbour -> Vertex 2
            %         11. South West Neighbour -> Vertex 3
            %         20. South Neighbour -> Vertex 4
            % Vertex 2
            %         20. South Neighbour -> Vertex 3
            %         21. South East Neighbour -> Vertex 4
            %         30. East Neighbour -> Vertex 1
            % Vertex 3
            %         30. East Neighbour -> Vertex 4
            %         31. North East Neighbour -> Vertex 1
            %         40. North Neighbour -> Vertex 2
            % Vertex 4
            %         40. North Neighbour -> Vertex 1
            %         41. North West Neighbour -> Vertex 2
            %         10. West Neighbour -> Vertex 3
            
            switch vertexNum
                case 1
                    Neighbours = [10,11,20];
                case 2
                    Neighbours = [20,21,30];
                case 3
                    Neighbours = [30,31,40];
                case 4
                    Neighbours = [40,41,10];
            end
            
            for ii = 1:3
                
                Neighbour = Neighbours(ii);
                dir = floor(Neighbour/10);
                Ntyp = mod(Neighbour,10);
                
                [exist_NQ, refNQ] = refNeighbour(refLeaf,dir,Ntyp);
                if exist_NQ == 1
                    
                    idxNQ = findNeighbour(Quadtree,idxLeaf,refLeaf, refNQ);
                    
                    refNQ = Quadtree.Node{idxNQ,1}{2,1}(1:end);
                    levelNQ = length(refNQ);
                    
                    if Quadtree.isleaf(idxNQ) && levelLeaf > levelNQ
                        
                        [Quadtree] = Decompose_helper(Quadtree,NURBS,idxNQ);
                    
                    end
                    
                end
            end
            
        end
    end
end


data = Quadtree.get(1);
oNURBS = data{3};

numLeaves = length(updatedLeaves);

for i = 1:numLeaves
    
    LEAF = updatedLeaves(i);
    
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
    
    data = Quadtree.Node{LEAF.idx,1};
    
    refLeaf = data{2};
    
    % seach Neighbours that share the same vertex
    % Vertex 1
    %         10. West Neighbour -> Vertex 2
    %         11. South West Neighbour -> Vertex 3
    %         20. South Neighbour -> Vertex 4
    % Vertex 2
    %         20. South Neighbour -> Vertex 3
    %         21. South East Neighbour -> Vertex 4
    %         30. East Neighbour -> Vertex 1
    % Vertex 3
    %         30. East Neighbour -> Vertex 4
    %         31. North East Neighbour -> Vertex 1
    %         40. North Neighbour -> Vertex 2
    % Vertex 4
    %         40. North Neighbour -> Vertex 1
    %         41. North West Neighbour -> Vertex 2
    %         10. West Neighbour -> Vertex 3
    
    switch LEAF.vertexNum
        case 1
            Neighbours = [10,11,20];
            updateVertices = [2,3,4];
            Children = [4,3,1];
        case 2
            Neighbours = [20,21,30];
            updateVertices = [3,4,1];
            Children = [3,1,2];
        case 3
            Neighbours = [30,31,40];
            updateVertices = [4,1,2];
            Children = [1,2,4];
        case 4
            Neighbours = [40,41,10];
            updateVertices = [1,2,3];
            Children = [2,4,3];
    end
    
    for ii = 1:3
        
        Neighbour = Neighbours(ii);
        updateVertex = updateVertices(ii);
        dir = floor(Neighbour/10);
        Ntyp = mod(Neighbour,10);
        
        [exist_NQ, refNQ] = refNeighbour(refLeaf,dir,Ntyp);
        if exist_NQ == 1
            
            idxNQ = findNeighbour(Quadtree,LEAF.idx,refLeaf, refNQ);
            
            if Quadtree.isleaf(idxNQ)
                idx = idxNQ;
            else
                idxNQchildren = Quadtree.getchildren(idxNQ);
                idx = idxNQchildren(Children(ii));
            end
            
            % get data stored in current leaf's Neighbour
            data = Quadtree.Node{idx,1};
            % intersection points
            intersectionPoints = data{3};
            % intersection points location
            locIntersectionPoints = data{4};
            % NURBS segment
            NURBS = data{6};
            % quad vertices 
            quadVertices = data{7};
            
            currentVertexCoor = quadVertices(:,updateVertex);
            
            if norm(currentVertexCoor - LEAF.oldVertexCoor(1:2)) < tol
                
                % replace quad vertex by its projection
                quadVertices(:,updateVertex) = LEAF.vertexCoor(1:2) ;
                
            end
            
            % check if current leaf's NB contains a NURBS segment
            if ~isempty(NURBS) % yes
                
                idxU1 = find( vecnorm(intersectionPoints-LEAF.oldIntersectionPoints(:,1) ) < tol);
                idxU2 = find( vecnorm(intersectionPoints-LEAF.oldIntersectionPoints(:,2) ) < tol);
                
                if ~isempty(idxU1)
                    idxU = idxU1;
                elseif ~isempty(idxU2)
                    idxU = idxU2;
                else
                    idxU = 0;
                end
                
                if idxU ~= 0
                    
                    % update intersection points
                    intersectionPoints(:,idxU) = LEAF.vertexCoor;
                    locIntersectionPoints(idxU) = updateVertex;
                    intersections = intersectionPoints(3,:);
                    intersections = sort(intersections);
                    
                    % sort intersectionPoints
                    if ( locIntersectionPoints(1) > locIntersectionPoints(2) )
                        locIntersectionPoints = locIntersectionPoints(end:-1:1);
                        intersectionPoints = intersectionPoints(:,end:-1:1);
                    end
                    
                    data{3} = intersectionPoints;
                    data{4} = locIntersectionPoints;
                    data{5} = intersections;
                    
                    newNURBS = oNURBS;
                    
                    for j = 1:length(intersections)
                        % loop over newKnotVals
                        % for knotVal in newknots
                        KnotVal = intersections(j);
                        numKnotIns = newNURBS.degree-sum(abs(newNURBS.knots(:) - KnotVal) < 1e-10);
                        % Insert knot numKnotsIns times
                        if numKnotIns > 0
                            [newNURBS.knots, newNURBS.controlPoints, newNURBS.weights] = CurveKnotIns(newNURBS.degree,...
                                newNURBS.knots, newNURBS.controlPoints, newNURBS.weights, KnotVal, numKnotIns);
                        end
                    end
                    
                    [newNURBS] = extractNURBS_segment(intersections, newNURBS);
                    
                    % update NURBS segment
                    data{6} = newNURBS;
                end
            else % no
                % update intersection points
                data{3} = LEAF.vertexCoor;
                data{5} = LEAF.vertexCoor(3);
            end
            
            % update quad vertices
            data{7} = quadVertices;
            % update data stored in current leaf's Neigbour
            Quadtree = Quadtree.set(idx, data);
            
            % plot quad (current leaf's Neighbour)
            if dbg == 1
                patch(quadVertices(1,:),quadVertices(2,:),'g','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
            end
        end
    end
end
end


