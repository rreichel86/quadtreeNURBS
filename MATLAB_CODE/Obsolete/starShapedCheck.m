function [K] = starShapedCheck(vertices)

%   starShapedCheck checks if a polygon is star-shaped
%   Returns the Kernel of the polygon K. If star-shaped K contains the
%   coordinates of the Kernel, otherwise K is empty

% Input:
% vertices: coordinates of the polygon to be searched for kernel
% Output:
% K: the coordinates of the kernel

% Check global orientation of polygon: -1=clockwise, 1=counter-clockwiese
oriP = starShapedCheckOrientation(vertices);

% Check convexity and concavity
% if oriP = ori(local) => convex

verticesClosed = horzcat(vertices(:,length(vertices)),vertices); % copy last vertex to the beginning of the array

% Initialize:
concaveVertices = [];
concaveVerticesCoords = [];
% Loop over the length of the vertices array in search for concave vertices
for i = 2:(length(verticesClosed)-1)
    vertex = i - 1;
    % Checks the orientation of each vertex against the global direction of
    % the polygon, if they don't matchm the vertex is concave
    verticesLocal = verticesClosed((1:2),(i - 1:i + 1));
    oriLocal = starShapedCheckOrientation(verticesLocal);
    if oriLocal ~= oriP
        concaveVertices = [concaveVertices, vertex];
        concaveVerticesCoords = [concaveVerticesCoords, (verticesLocal((1:2),(2:2)))];
        
    end
    
end

% Initialize K, at the end K will contain only the vertices of the kernel
K = vertices;

% If there are no concave vertices, the polygon coincides with its own kernel
if isempty(concaveVertices)
    return
end

KFirsttoEnd = [K,K(:,1)]; % Copy first vertex to the end of the array

%insert points
for v = 1:(length(concaveVertices))
    % V = coordinates of current vertex
    V = concaveVerticesCoords(:,v);
    % VBack = coordinates of previous vertex
    VBack = verticesClosed(:,(concaveVertices(:,v)));
    VForwardMatrix = [verticesClosed, verticesClosed(:,2)];
    % VFroward = coordinates of next vertex
    VForward = VForwardMatrix(:,(concaveVertices(:,v)+2));
    insertions = 0;
    deletion = 0;
    s = 1;
    % for the elongation of each of the convex vertices, check if there
    % is an intersection point with one of the sides of the polygon,
    % according to the algorithm proposed by Zhao et al (2010),
    % DOI: 10.1109/ICCASM.2010.5620799
    while s <= (length(K))
        S = K(:,s);
        SForward = KFirsttoEnd(:,s+1);
        tri1 = [VBack,V,S];
        tri2 = [VBack,V,SForward];
        tri3 = [V,VForward,S];
        tri4 = [V,VForward,SForward];
        oriTri1 = starShapedCheckOrientation(tri1);
        oriTri2 = starShapedCheckOrientation(tri2);
        oriTri3 = starShapedCheckOrientation(tri3);
        oriTri4 = starShapedCheckOrientation(tri4);
        
        if oriTri1 == -1*oriTri2 && oriTri1 ~= 0 && oriTri2 ~= 0
            L1 = [VBack,V];
            L2 = [S,SForward];
            u = starShapedIntersect(L1,L2); % obtain the coordinats of the intersection
            K = [K(:,1:s), u, K(:,s+1:end)]; % add them to K
            KFirsttoEnd = [K,K(:,1)];
            insertions = insertions + 1;
        elseif oriTri3 == -1*oriTri4 && oriTri3 ~= 0 && oriTri4 ~= 0
            L1 = [V,VForward];
            L2 = [S,SForward];
            u2 = starShapedIntersect(L1,L2); % obtain the coordinats of the intersection
            K = [K(:,1:s), u2, K(:,s+1:end)]; % add them to K
            KFirsttoEnd = [K,K(:,1)];
            insertions = insertions + 1;
        end
        s = s + 1;
    end
    
    % Delete the vertices outside of the resulting kernel
    for s = 1:(length(K))
        % if K has more than 2 point, hence isn't a line
        if length(K) > 2
            VBack = verticesClosed(:,(concaveVertices(:,v)));
            VForwardMatrix = [verticesClosed, verticesClosed(:,2)];
            VForward = VForwardMatrix(:,(concaveVertices(:,v)+2));
            S = K(:,s-deletion);
            tri1 = [VBack,V,S];
            tri3 = [V,VForward,S];
            oriTri1 = starShapedCheckOrientation(tri1);
            oriTri3 = starShapedCheckOrientation(tri3);
            if oriP == 1 && (oriTri1 == -1 || oriTri3 == -1)
                K = [K(:,1:s-deletion-1), K(:,s-deletion+1:end)];
                KFirsttoEnd = [K,K(:,1)];
                deletion = deletion + 1;
            end
            if oriP == -1 && (oriTri1 == 1 || oriTri3 == 1)
                K = [K(:,1:s-deletion-1), K(:,s-deletion+1:end)];
                KFirsttoEnd = [K,K(:,1)];
                deletion = deletion + 1;
                
            end
            
        else
            K = [];
        end
    end
    
    
    
end

% Check if the kernel isn't a line or has too small area
if length(K) < 3 || round(polyarea(K(1,:), K(2,:)), 9) == 0
    K = [];
    return
end

end

