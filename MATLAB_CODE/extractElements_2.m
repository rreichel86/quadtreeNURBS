function [numcoor,coor,numel,connectivity,maxnel,...
          numKnotVectors,knotVectors,idxControlPoints] = extractElements_2(Quadtree)
% extractElements: get polygonal elements from Quadtree data structure
%
% INPUT:
% Quadtree Data
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

%% 
% Get NURBS curve 
data = Quadtree.Node{1,1};
NURBS.degree = data{3};
NURBS.knots  = data{4};
NURBS.controlPoints = data{5};
NURBS.weights = data{6};

% number of knots
NURBS_nknots = length(NURBS.knots);

% number of control points
NURBS_ncp = length(NURBS.controlPoints);



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

j = 1;
% counter for the number of knot vectors
countKnotVectors = 0;
% loop over leaves 
for i = 1:numleaves
    
     % current Quad reference
    refLeaf = Quadtree.Node{leaves(i),1}{2,1}(1:end);
    intersections = Quadtree.Node{leaves(i),1}{5,1};
    
    % Check neighborhood of current leaf
    % get mid points if they exist
    [midPoints] = getMidPoints(Quadtree,leaves(i),refLeaf);

    if isempty(intersections) || length(intersections) == 1
        
        % quad's geometry
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        % intersection points
        intersectionPoints = [];

        
        coordinates = [quad;...       % nodal x-coor and y-coor
                       ones(1,4);...  % weights
                       ones(1,4);...  % type 
                       zeros(1,4);... % which_region
                       -1*ones(1,4)]; % inside_region
                   
        element = [quad,midPoints];
        
    elseif isempty(Quadtree.Node{leaves(i),1}{3,1}) == 1
        
        % quad's geometry 
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        % intersection points
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1), Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        
        % counter for the number of knot vectors
        countKnotVectors = countKnotVectors + 1;
        % intersection points weights 
        weights = [Quadtree.Node{leaves(i),1}{9,1}(:,1), Quadtree.Node{leaves(i),1}{9,1}(:,end)];
        
        knotVectors{countKnotVectors} = [countKnotVectors,NURBS.degree,intersections,NURBS_nknots,NURBS.knots];
        controlPoints_coor{countKnotVectors} = [NURBS.controlPoints',NURBS.weights'];
        
        idxControlPoints{countKnotVectors} = zeros(1,NURBS_ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = NURBS_ncp;

        coordinates = [     quad, intersectionPoints, NURBS.controlPoints;... % nodal/ctrlPts x-coor and y-coor
                        ones(1,4),           weights,       NURBS.weights;... % weights
                        ones(1,4),       2*ones(1,2), 2*ones(1,NURBS_ncp);... % type
                        zeros(1,4),       zeros(1,2),  zeros(1,NURBS_ncp);... % which_region
                       -1*ones(1,4),      zeros(1,2),  zeros(1,NURBS_ncp)];   % inside_region
                   
        element = [quad,midPoints,intersectionPoints];
        
    elseif isempty(Quadtree.Node{leaves(i),1}{4,1}) == 1
        
        % quad's geometry 
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        % intersection points
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1), Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        
        % intersection points weights 
        weights = [Quadtree.Node{leaves(i),1}{9,1}(:,1), Quadtree.Node{leaves(i),1}{9,1}(:,end)];
        
        % counter for the number of knot vectors
        countKnotVectors = countKnotVectors + 1;
        
        knotVectors{countKnotVectors} = [countKnotVectors,NURBS.degree,intersections,NURBS_nknots,NURBS.knots];
        controlPoints_coor{countKnotVectors} = [NURBS.controlPoints',NURBS.weights'];
        
        idxControlPoints{countKnotVectors} = zeros(1,NURBS_ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = NURBS_ncp;

        coordinates = [     quad, intersectionPoints, NURBS.controlPoints;... % nodal/ctrlPts x-coor and y-coor
                        ones(1,4),           weights,       NURBS.weights;... % weights
                        ones(1,4),       2*ones(1,2), 2*ones(1,NURBS_ncp);... % type
                        zeros(1,4),       zeros(1,2),  zeros(1,NURBS_ncp);... % which_region
                       -1*ones(1,4),      zeros(1,2),  zeros(1,NURBS_ncp)];   % inside_region
                   
        element = [quad,midPoints,intersectionPoints];
        
    else

        % quad's geometry 
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        % intersection points
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1), Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        
        % intersection points weights 
        weights = [Quadtree.Node{leaves(i),1}{9,1}(:,1), Quadtree.Node{leaves(i),1}{9,1}(:,end)];
        
        % counter for the number of knot vectors
        countKnotVectors = countKnotVectors + 1;
        
        knotVectors{countKnotVectors} = [countKnotVectors,NURBS.degree,intersections,NURBS_nknots,NURBS.knots];
        controlPoints_coor{countKnotVectors} = [NURBS.controlPoints',NURBS.weights'];
        
        idxControlPoints{countKnotVectors} = zeros(1,NURBS_ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = NURBS_ncp;

        coordinates = [     quad, intersectionPoints, NURBS.controlPoints;... % nodal/ctrlPts x-coor and y-coor
                        ones(1,4),           weights,       NURBS.weights;... % weights
                        ones(1,4),       2*ones(1,2), 2*ones(1,NURBS_ncp);... % type
                        zeros(1,4),       zeros(1,2),  zeros(1,NURBS_ncp);... % which_region
                       -1*ones(1,4),      zeros(1,2),  zeros(1,NURBS_ncp)];   % inside_region
                   
        element = [quad,midPoints,intersectionPoints];
        
    end
    
    intersections_coor{i} = intersectionPoints';% intersection points
    coor{i} = coordinates;% nodal coordinates
    
    % remove -99 values and arrange the element nodal coordinates
    % counterclockwise starting from the left bottom corner
    element = element(element ~= -99);
    ncol = size(element, 1);
    element = reshape(element,[2,ncol/2]);
    x = element(1,:);
    y = element(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);% angle counterclockwise from [-pi,pi]
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order);
    if a(1,1) > min(a(1,2:end))
        %if there is any point angle which is less than the left bottom
        %point than it should be end point
        element(1,:) = [x_1(1,2:end),x_1(1,1)];
        element(2,:) = [y_1(1,2:end),y_1(1,1)];
    else
        element(1,:) = x_1;
        element(2,:) = y_1;
    end
    %if int_cor empty than there is only one element that was itself leaf
    %otherwise it would be two elements
    if isempty(intersectionPoints)    
        elements{j} = element;
        % element number
        connectivity{j}(1,1) = j;
        % knot vector number
        connectivity{j}(1,2) = 0;
        % region number 
        connectivity{j}(1,3) = 1;
        j = j+1;
    else
        % loop over intersection points 
        for k = 1:2
            % search for intersection points in current element nodal coordinates
            [b] = find ( abs(intersectionPoints(1,k) - element(1,:))<1e-10);
            [d] = find (abs(intersectionPoints(2,k) - element(2,:))<1e-10);
            col(:,k) = intersect(b,d);
        end
        % for the control points in Clockwise direction 
        if col(1,1) < col(1,2)
            elements{j} = [element(:,[1:col(1,1)]),element(:,[col(1,2):end])];
            elements{j+1} = [element(:,[col(1,1):col(1,2)])];
        else
            elements{j} = [element(:,[1:col(1,2)]),element(:,[col(1,1):end])];
            elements{j+1} = [element(:,[col(1,2):col(1,1)])];   
        end
        % element numbers
        connectivity{j}(1,1) = j;
        connectivity{j+1}(1,1) = j+1;
        % knot vector number
        connectivity{j}(1,2) = countKnotVectors;
        connectivity{j+1}(1,2) = countKnotVectors;
        % region number 
        connectivity{j}(1,3) = 1;
        connectivity{j+1}(1,3) = 1;
        j = j+2;
    end
end % end loop over leaves

%%
% intersection point coordinates in matrix form
intersections_coor = cell2mat(intersections_coor);
tol = 1e-10;

% coordinates in matrix form
coor = [coor{:}]';
% remove repeated coordinates and rearrange them by increasing x-coor value
tmp_coor = uniquetol(coor,tol,'ByRows',true);

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
    kernel = computePolygonKernel( elements{iel}' );
    % compute scaling center
    % compute centroid of the kernel
    [scx,scy] = centroid(polyshape(kernel,'Simplify',false));
    coor(numcoor0+iel,2:7) = [scx,scy,1,1,0,-1];
end

% coor numbers 
coor(:,1) = (1:numcoor);

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

indices = zeros(size(intersections_coor,1),1);
% loop over intersection points
for k = 1:size(intersections_coor,1)
    % search for intersection point in coor 
    % and get its index
    [a] = find ( abs(coor(:,2)-intersections_coor(k,1))<1e-10);
    [b] = find ( abs(coor(:,3)-intersections_coor(k,2))<1e-10);
    index = intersect(a,b);
    indices(k,1) = index;
end

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

numKnotVectors = countKnotVectors;
numIdxControlPoints = countKnotVectors;
for ictrlp = 1:numIdxControlPoints
    ctrlp = controlPoints_coor{ictrlp};
    nctrlp = idxControlPoints{ictrlp}(1,2);
    
    for n = 1 : nctrlp
        
        a = find ( abs(coor(:,2)-ctrlp(n,1))<1e-10);
        b = find ( abs(coor(:,3)-ctrlp(n,2))<1e-10);
        index = intersect(a,b);
        idxControlPoints{ictrlp}(1,n+2) = index';
    end   
end

end




