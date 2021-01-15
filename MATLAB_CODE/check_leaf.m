function [Quadtree] = check_leaf(Quadtree)
% plot_leaf: plots each quadtree's leaf


% figure(2)
% axis square
% hold on

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
    
    % Intersections
    intersections = data{5};
    horzIntrsc = data{3}';
    vertIntrsc = data{4}';
    
    intrsc = [horzIntrsc;vertIntrsc];
    
    % check if leaf has intersections
    if isempty(intersections) || length(intersections) == 1
        continue
    else
        
        % Quad definition
        Quad = data{7};
        
        Qxmin = min(Quad(1,:));
        Qymin = min(Quad(2,:));
        Qxmax = max(Quad(1,:));
        Qymax = max(Quad(2,:));
        
        if ~isempty(horzIntrsc) && ~isempty(vertIntrsc)
            
            eta_h = (intrsc(1,1) - Qxmin)/( Qxmax - Qxmin );
            eta_v = (intrsc(2,2) - Qymin)/( Qymax - Qymin );
            
            % horizontal intersection with quad's
            if abs(intrsc(1,2) - Qymin) < tol %  bottom edge
                vertices = [1,2];
            else % top edge
                vertices = [4,3];
                eta_v = 1-eta_v;
            end
            
            % vertical intersection with quad's
            if abs(intrsc(2,1) - Qxmin) < tol % left edge
                vertexNum = vertices(1);
            else % right edge
                vertexNum = vertices(2);
                eta_h = 1-eta_h;
            end
            
            if eta_h < 1/5 && eta_v < 1/5
                
                n = [-(intrsc(2,2)-intrsc(1,2));...
                    intrsc(2,1)-intrsc(1,1)];
                
                n = n/norm(n);
                p = intrsc(1,:)';
                v = Quad(:,vertexNum);
                newVertexCoor = v + n.'*(p-v)*n;
                
                oldVertexCoor = Quad(:,vertexNum);
                Quad(:,vertexNum) = newVertexCoor;
                
%                 % plot vertex
%                 plot(Quad(1,vertexNum),Quad(2,vertexNum), 'o','Color','red','MarkerFaceColor','r','MarkerSize',6);
%                 % plot quad 
%                 patch(Quad(1,:),Quad(2,:),'r','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
                
                % update quad definition
                data{7} = Quad;
                data{5} = [];
                Quadtree = Quadtree.set(idxLeaf, data);
                
                LEAF.idx = idxLeaf;
                LEAF.vertexNum = vertexNum;
                LEAF.vertexCoor = newVertexCoor;
                LEAF.oldVertexCoor = oldVertexCoor;
                
                updatedLeaves = [updatedLeaves LEAF];
                
            end
        else
            continue
        end
    end
end

numLeaves = length(updatedLeaves);

for i = 1:numLeaves
    
    LEAF = updatedLeaves(i);
%     idxLeaf = LEAF.idx;
    
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
    
    % Quad definition
%     Quad = data{7};

    % plot vertex
%     plot(Quad(1,LEAF.vertexNum),Quad(2,LEAF.vertexNum), 'o','Color','red','MarkerFaceColor','r','MarkerSize',6);
    % plot quad
%     patch(Quad(1,:),Quad(2,:),'r','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
    
%     Qxmin = min(Quad(1,:));
%     Qymin = min(Quad(2,:));
%     Qxmax = max(Quad(1,:));
%     Qymax = max(Quad(2,:));
    
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
            
            % quad definition 
            data = Quadtree.Node{idx,1};
            Quad = data{7};
            
%             % plot quad 
%             patch(Quad(1,:),Quad(2,:),'g','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
            
            currentVertexCoor = Quad(:,updateVertex);
            
            if norm(currentVertexCoor - LEAF.oldVertexCoor) < tol
                
                Quad(:,updateVertex) = LEAF.vertexCoor;
                % update quad definition
                data{7} = Quad;
                Quadtree = Quadtree.set(idx, data);
            end 
            
        end
    end
end
end


