function [numcoor,coor,numel,connectivity,maxnel,...
          numKnotVectors,knotVectors,idxControlPoints] = extractElements_2(Quadtree)
% extractElements: get polygonal elements from Quadtree data structure
%
% INPUT:
% Quadtree ------------------- Quadtree data structure
%
% OUTPUT:
% numcoor -------------------- number of coordinates = number of nodes
% coor ----------------------- nodes coordinates and weights 
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary
%                                             1 - inside 
%                                            -1 - outside
%
% numel ---------------------- number of elements
% connectivity --------------- elements connectivity matrix as nel-tupel of 
%                              nodes, where the first six entries
%                              iel - element number
%                              ikv - knot vector number 
%                              which_region - region number
%                              inode - index initial node 
%                              jnode - index last node
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, which_region, inode, jnode, nel, node_1,...,node_nel, scaling_center]
% maxnel --------------------- maximum number of nodes on any element
%
% numKnotVectors ------------- number of knot vectors
% knotVectors ---------------- contains knot vectors and following information
%                              ikv - knot vektor number
%                              degree - NURBS curve degree
%                              iknot - initial knot value
%                              jknot - end knot value
%                              nkonts - number of knots per knot vector
%
% knotVectors = [ikv, degree, iknot, jknot, nknots, knot_1,...,knot_nknots]
% idxControlPoints ----------- control points indices
% idxControlPoints = [icp, ncp, idx_1,...idx_ncp] 
%
% -----------------------------------------------------------------------------

tol = 1e-10;
%% 
% Get NURBS curve 
data = Quadtree.Node{1,1};
NURBS = data{3};

% compute point of the NURBS curve
NURBS_pts = CalculateNURBS(NURBS);

% Compute bounding box that enclosed the NURBS curve
x_min = min(NURBS_pts(:,1));
x_max = max(NURBS_pts(:,1));
y_min = min(NURBS_pts(:,2));
y_max = max(NURBS_pts(:,2));

%%
% Get Quadtree leaves
leaves = Quadtree.findleaves();
numleaves = length(leaves);

% Get total number of elements
[numel] = countElements(Quadtree,leaves);

%%
% prealloc cell array
coor = cell(numleaves,1);
controlPoints_coor = cell(numleaves,1);
knotVectors = cell(numleaves,1);
elements = cell(numel,1);
connectivity = cell(numel,1);
intersections_coor = cell(numleaves,1);
idxControlPoints = cell(numleaves,1);

% counter for the number of elements
ielno = 1;
% counter for the number of knot vectors
countKnotVectors = 0;
% loop over leaves 
for i = 1:numleaves
    
    idxLeaf = leaves(i);
    
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
    
    % leaf reference
    refLeaf = data{2};
    % quad vertices
    quadVertices = data{7};
    % number of quad vertices
    numQuadVertices = size(quadVertices,2);
    % intersecion points
    intersectionPoints = data{3};
    % number of intersection points
    numIntersectionPoints = size(intersectionPoints,2);
    % intersection points location
    locIntersectionPoints = data{4};
    
    % Check neighborhood of current leaf
    % get mid points if they exist
    [numMidPoints,midPoints,locMidPoints] = getMidPoints(Quadtree,idxLeaf,refLeaf);
    
    % number of polygon's vertices
    numPolyVertices = numQuadVertices + numMidPoints + numIntersectionPoints;
    % polygon vertices
    poly = zeros(2, numPolyVertices);
    % intersection points indices
    idxIntersectionPoints = zeros(1,numIntersectionPoints);
    
    % arrange quadVertices, midPoints and intersectionPoints CCW in poly
    iMp = 0;
    iIp = 0;
    nPv = 0;
    for iQv = 1:numQuadVertices
        % current edge's endpoint
        A = quadVertices(:,iQv);
        % insert quad's edge endpoint
        nPv = nPv + 1;
        poly(:,nPv) = A;
        % check if a mid point lies on the current edge
        if numMidPoints ~= 0
            iMp = find(locMidPoints == iQv);
            if isempty(iMp)
                iMp = 0;
            end
        end
        % check if an intersection point lies on the current edge
        if numIntersectionPoints ~= 0
            iIp = find(locIntersectionPoints == iQv);
            if isempty(iIp)
                iIp = 0;
            end
        end
        % a mid point and an intersection point lie on the current edge
        if iMp ~= 0 && iIp ~= 0
            
            % mid point
            M = midPoints(:,iMp);
            % intersection point
            P = intersectionPoints(1:2,iIp);
            
            % distance to current edge's endpoint
            dist_AP = norm(P-A);
            dist_AM = norm(M-A);
            dist_MP = norm(P-M);
            
            % is P closer to A than M?
            if dist_AP < dist_AM
                % check if P == A?
                if dist_AP < tol
                    idxIntersectionPoints(iIp) = nPv;
                else
                    % insert intersection point
                    nPv = nPv + 1;
                    poly(:,nPv) = P;
                    idxIntersectionPoints(iIp) = nPv;
                end
                % insert mid point
                nPv = nPv + 1;
                poly(:,nPv) = M;
            else
                % insert mid point
                nPv = nPv + 1;
                poly(:,nPv) = M;
                % check if P == M ?
                if dist_MP < tol
                    idxIntersectionPoints(iIp) = nPv;
                else
                    % insert intersection point
                    nPv = nPv + 1;
                    poly(:,nPv) = P;
                    idxIntersectionPoints(iIp) = nPv;
                end
            end
            
            %  a mid point lies on the current edge
        elseif  iIp ~= 0
            
            % intersection point
            P = intersectionPoints(1:2,iIp);
            
            % distance to current edge's endpoint
            dist_AP = norm(P-A);
            
            % check if P == A?
            if dist_AP < tol
                idxIntersectionPoints(iIp) = nPv;
            else
                % insert intersection point
                nPv = nPv + 1;
                poly(:,nPv) = P;
                idxIntersectionPoints(iIp) = nPv;
            end
            
            % a mid point lies on the current edge
        elseif  iMp ~= 0
            
            % mid point
            M = midPoints(:,iMp);
            % insert mid point
            nPv = nPv + 1;
            poly(:,nPv) = M;
        end
    end
    
    % update number of polygon's vertices
    if (nPv ~= numPolyVertices)
        numPolyVertices = nPv;
    end
    
    % check if leaf has intersections
    if numIntersectionPoints == 0
        
        % intersection points
        intersectionPoints = [];
        % nodal coordinates
        coordinates = zeros(6, 4);
         % quad vertices
        coordinates(1:2,1:4) = quadVertices; % nodal x-coor and y-coor
        coordinates(4,1:4) = 1;       % type  
        
%     elseif  numIntersectionPoints == 1
%          
%         % intersection points
%         intersectionPoints = [];
%         % nodal coordinates
%         coordinates = zeros(6, 4);
%          % quad vertices
%         coordinates(1:2,1:4) = quadVertices; % nodal x-coor and y-coor
%         coordinates(4,1:4) = 1;       % type 
    else
        % degree
        degree = NURBS.degree;
        % knot vector
        knots = NURBS.knots;
        % control points
        controlPoints = NURBS.controlPoints;
        % weights
        weights = NURBS.weights;
        
        if numIntersectionPoints == 1
            % intersections
            iknot = data{5}(1);
            jknot = data{5}(1);
            % intersection points
            intersectionPoints = [intersectionPoints(1:2,:),intersectionPoints(1:2,:)];
            
        elseif numIntersectionPoints == 2
            % intersections
            iknot = data{5}(1);
            jknot = data{5}(2);
            % intersection points
            intersectionPoints = intersectionPoints(1:2,:);
            
        end
        % number of knots
        nknots = length(knots);
        % number of control points / weights
        ncp = length(controlPoints);
        
        % counter for the number of knot vectors
        countKnotVectors = countKnotVectors + 1;
        knotVectors{countKnotVectors} = [countKnotVectors,degree,iknot,jknot,nknots,knots];
        controlPoints_coor{countKnotVectors} = [controlPoints',weights'];
        
        idxControlPoints{countKnotVectors} = zeros(1,ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = ncp;
        
        % nodal coordinates
        coordinates = zeros(6, 6 + ncp);
        % quad vertices 
        coordinates(1:2,1:4) = quadVertices; % nodal x-coor and y-coor
        coordinates(4,1:4) = 1;       % type
        % intersection points 
        coordinates(1:2,5:6) = intersectionPoints; % nodal x-coor and y-coor
        coordinates(4,5:6) = 2;       % type
        % control points 
        coordinates(1:2,7:6+ncp) = controlPoints; % ctrlPts x-coor and y-coor
        coordinates(4,7:6+ncp)   = 2; % type                   
    end
    
    intersections_coor{i} = intersectionPoints';% intersection points
    coor{i} = coordinates;% nodal coordinates
    
    % check if leaf has intersections
    % leaf is splitted in 2 elements if it has 2 intersection points
    if numIntersectionPoints == 0 
        elements{ielno} = poly(1:2,1:numPolyVertices);
        % element number
        connectivity{ielno}(1,1) = ielno;
        % knot vector number
        connectivity{ielno}(1,2) = 0;
        % region number
        connectivity{ielno}(1,3) = 1;
        ielno = ielno+1;
    elseif  numIntersectionPoints == 1
        elements{ielno} = poly(1:2,1:numPolyVertices);
        % element number
        connectivity{ielno}(1,1) = ielno;
        % knot vector number
        connectivity{ielno}(1,2) = countKnotVectors;
        % region number
        connectivity{ielno}(1,3) = 1;
        ielno = ielno+1;
    else
        % split leaf in 2 elements 
        elements{ielno} = [poly(1:2,1:idxIntersectionPoints(1)),...
            poly(1:2,idxIntersectionPoints(2):numPolyVertices)];
        
        elements{ielno+1} = [poly(1:2,idxIntersectionPoints(1):idxIntersectionPoints(2))];
        
        % element numbers
        connectivity{ielno}(1,1) = ielno;
        connectivity{ielno+1}(1,1) = ielno+1;
        % knot vector number
        connectivity{ielno}(1,2) = countKnotVectors;
        connectivity{ielno+1}(1,2) = countKnotVectors;
        % region number
        connectivity{ielno}(1,3) = 1;
        connectivity{ielno+1}(1,3) = 1;
        ielno = ielno+2;
    end
end % end loop over leaves

%%
% intersection point coordinates in matrix form
intersections_coor = cell2mat(intersections_coor);

% coordinates in matrix form
coor = [coor{:}]';
% remove repeated coordinates and rearrange them by increasing x-coor value
[~, idxTemp, ~] = uniquetol(coor(:,1:2),tol,'ByRows',true);
tmp_coor = coor(idxTemp,:);

% number of coordinates without counting scaling centers
numcoor0 = size(tmp_coor,1);
% number of scaling centers = number of elements
numsc = numel;
% update numcoor
numcoor = numcoor0 + numsc;

% prealloc coor matrix 
coor = zeros(numcoor, 7);
% coordinates 
coor(1:numcoor0,2:7) = tmp_coor;

% compute scaling center of polygonal elements
for iel = 1:numel
    % compute kernel
    kernel = computePolygonKernel( elements{iel}' ,0);
    % compute scaling center
    % compute centroid of the kernel
    [scx,scy] = centroid(polyshape(kernel,'Simplify',false));
    coor(numcoor0+iel,2:7) = [scx,scy,1,1,0,-1];
end

% coor numbers 
coor(:,1) = (1:numcoor);

numKnotVectors = countKnotVectors;
numIdxControlPoints = countKnotVectors;
% loop over control points
for ictrlp = 1:numIdxControlPoints
    ctrlp = controlPoints_coor{ictrlp};
    nctrlp = idxControlPoints{ictrlp}(1,2);
    
    for n = 1 : nctrlp
        % search for control points in coor 
        % and get its index
        index = find ( vecnorm(coor(:,2:3)-ctrlp(n,1:2),2,2)<tol);
        idxControlPoints{ictrlp}(1,n+2) = index;
        
        % update nodes weights and type
        coor(index,4) = ctrlp(n,3);
        coor(index,5) = 2;
    end   
end

indices = zeros(size(intersections_coor,1),1);
% loop over intersection points
for k = 1:size(intersections_coor,1)
    % search for intersection point in coor
    % and get its index
    index = find( vecnorm(coor(:,2:3)-intersections_coor(k,1:2),2,2) < tol);
    indices(k,1) = index;
    
    % update type
    coor(index,5) = 2;
end

coor(coor(:,5) == 1,7) = -1;

% loop over nodes, excluding control points
for ii = find(coor(:,5) == 1)'
    % Check if current node is inside the bounding box
    if ( isPointInQuad([x_min,y_min], [x_max,y_max], coor(ii,2:3)) == 1 )
        % Check if current node is also inside 
        % the region enclosed by the NURBS curve
        pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor(ii,2:3));
        coor(ii,7) = pointInPoly;
    end
end

% delete tmp_coor
clearvars tmp_coor



maxnel = 0;
% loop over elements
for iel = 1:numel
    % number of nodes per element
    nel = size(elements{iel},2); 
    % search for element nodal coordinates in coor
    % and get the corresponding node numbers
    for n = 1 : nel
        [a] = find ( abs(coor(:,2) - elements{iel}(1,n))<1e-10);
        [b] = find ( abs(coor(:,3) - elements{iel}(2,n))<1e-10);
        nodeNumber = intersect(a,b);   
        connectivity{iel}(1,6+n) = nodeNumber;
    end
    % remove repeated values without changing the order
    connectivity{iel}(1,7:end) = unique(connectivity{iel}(1,7:end),'stable');
    % update number of nodes per element
    nel = size(connectivity{iel}(7:end),2); 
    
    if sum(coor( connectivity{iel}(1,7:6+nel), 7)) < 0 % outside
        connectivity{iel}(1,3) = 1;
    else % inside
        connectivity{iel}(1,3) = 2;
    end
    
    connectivity{iel}(1,6) = nel;
    % maximum number of nodes on any element
    maxnel = max(nel,maxnel);  
    % append correspondig scaling center index
    connectivity{iel}(1,end+1) = numcoor0+iel;  
    
    ikv = connectivity{iel}(1,2);
    
    if ikv ~= 0
        
        elmt = connectivity{iel}(1,7:end);
        ecoor = coor( elmt(1:end), 2:3);
        inode = find( elmt == indices(2*ikv-1) );
        jnode = find( elmt == indices(2*ikv) );
       
        OP = orientation( ecoor(end,:), ecoor(inode,:) ,ecoor(jnode,:) );
        
        if OP == 1
            connectivity{iel}(1,4) = inode;
            connectivity{iel}(1,5) = jnode;
        else     
            connectivity{iel}(1,4) = jnode;
            connectivity{iel}(1,5) = inode;
        end 
    end 
end

end




