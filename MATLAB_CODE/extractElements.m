function [numcoor,coor,numel,connectivity,maxnel,...
          numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree)
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
%                              nodes, where the first three entries
%                              iel - element number
%                              ikv - knot vector number
%                              which_region - region number
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, which_region, nel, node_1,...,node_nel, scaling_center]
% maxnel --------------------- maximum number of nodes on any element
%
% numKnotVectors ------------- number of knot vectors
% knotVectors ---------------- contains knot vectors and following information
%                              ikv - knot vektor number
%                              degree - NURBS curve degree
%                              icp - index initial control point
%                              ecp - last control point
%                              nkonts - number of knots per knot vector
%
% knotVectors = [ikv, degree, icp, ecp, nknots, knot_1,...,knot_nknots]
% maxnknots ------------------ maximun number of knot on any knot vector
% idxControlPoints ----------- control points indices
% idxControlPoints = [icp, ncp, idx_1,...idx_ncp] 
%
% -----------------------------------------------------------------------------

references = cellfun(@(Q) Q(2),Quadtree.Node);
references{1} = [];

% Get Quadtree leaves
leaves = Quadtree.findleaves();
numleaves = length(leaves);

% Get total number of elements
[numel] = countElements(Quadtree,leaves);

% prealloc cell array
coor = cell(numleaves,1);
controlPoints_coor = cell(numleaves,1);
knotVectors = cell(numleaves,1);
elements = cell(numel,1);
connectivity = cell(numel,1);
intersections_coor = cell(numleaves,1);
idxControlPoints = cell(numleaves,1);

j = 1;
countKnotVectors = 0;
maxnknots = 0;% maximun number of knot on any knot vector
% loop over leaves 
for i = 1:numleaves
    
     % current Quad reference
    refLeaf = Quadtree.Node{leaves(i),1}{2,1}(1:end);
    intersections = Quadtree.Node{leaves(i),1}{5,1};
    
    % Check neighborhood of current leaf
    % get mid points if they exist
    [midPoints] = getMidPoints(Quadtree,references,refLeaf);

    if isempty(intersections) || length(intersections) == 1
        
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        % cp = [];% control points and weight
        intersectionPoints = [];

        coordinates = [quad;...       % nodal x-coor and y-coor
                       ones(1,4);...  % weights
                       ones(1,4)];    % type 
        element = [quad,midPoints];
        
    elseif isempty(Quadtree.Node{leaves(i),1}{3,1}) == 1
        
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        controlPoints = Quadtree.Node{leaves(i),1}{7,1};
        weights = Quadtree.Node{leaves(i),1}{9,1};
        ncp = length(controlPoints);
        % cp = [controlPoints',weights'];% control points and weights
        
        countKnotVectors = countKnotVectors + 1;
        knots = Quadtree.Node{leaves(i),1}{8,1};
        nknots = length(knots);
        maxnknots = max(maxnknots,nknots);
        degree = Quadtree.Node{leaves(i),1}{6,1};
        knotVectors{countKnotVectors} = [countKnotVectors,degree,0,0,nknots,knots];
        controlPoints_coor{countKnotVectors} = [controlPoints',weights'];
        idxControlPoints{countKnotVectors} = zeros(1,ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = ncp;
        
        coordinates = [quad, controlPoints;...   % nodal/ctrlPts x-coor and y-coor
                       ones(1,4), weights;...    % weights
                       ones(1,4), 2*ones(1,ncp)];  % type
                   
        element = [quad,midPoints,intersectionPoints];
        
    elseif isempty(Quadtree.Node{leaves(i),1}{4,1}) == 1
        
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        controlPoints = Quadtree.Node{leaves(i),1}{7,1};
        weights = Quadtree.Node{leaves(i),1}{9,1};
        ncp = length(controlPoints);
        % cp = [controlPoints',weights'];% control points and weights
        
        countKnotVectors = countKnotVectors + 1;
        knots = Quadtree.Node{leaves(i),1}{8,1};
        nknots = length(knots);
        maxnknots = max(maxnknots,nknots);
        degree = Quadtree.Node{leaves(i),1}{6,1};
        knotVectors{countKnotVectors} = [countKnotVectors,degree,0,0,nknots,knots];
        controlPoints_coor{countKnotVectors} = [controlPoints',weights'];
        idxControlPoints{countKnotVectors} = zeros(1,ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = ncp;

        coordinates = [quad, controlPoints;...   % nodal/ctrlPts x-coor and y-coor
                       ones(1,4), weights;...    % weights
                       ones(1,4), 2*ones(1,ncp)];  % type
                   
        element = [quad,midPoints,intersectionPoints];
        
    else

        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        
        controlPoints = Quadtree.Node{leaves(i),1}{7,1};
        weights = Quadtree.Node{leaves(i),1}{9,1};
        ncp = length(controlPoints);
        % cp = [controlPoints',weights'];% control points and weights
        
        countKnotVectors = countKnotVectors + 1;
        knots = Quadtree.Node{leaves(i),1}{8,1};
        nknots = length(knots);
        maxnknots = max(maxnknots,nknots);
        degree = Quadtree.Node{leaves(i),1}{6,1};
        knotVectors{countKnotVectors} = [countKnotVectors,degree,0,0,nknots,knots];
        controlPoints_coor{countKnotVectors} = [controlPoints',weights'];
        idxControlPoints{countKnotVectors} = zeros(1,ncp+2);
        idxControlPoints{countKnotVectors}(1,1) = countKnotVectors;
        idxControlPoints{countKnotVectors}(1,2) = ncp;
        
        coordinates = [quad, controlPoints;...   % nodal/ctrlPts x-coor and y-coor
                       ones(1,4), weights;...    % weights
                       ones(1,4), 2*ones(1,ncp)];  % type
                   
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
            elements{j} = [element(:,[1:col(1,1)]),controlPoints(:,[2:end-1]),element(:,[col(1,2):end])];
            elements{j+1} = [element(:,[col(1,1):col(1,2)]),fliplr(controlPoints(:,[2:end-1]))];
            
            
        else
            elements{j} = [element(:,[1:col(1,2)]),fliplr(controlPoints(:,[2:end-1])),element(:,[col(1,1):end])];
            elements{j+1} = [element(:,[col(1,2):col(1,1)]),controlPoints(:,[2:end-1])];   
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
coor(1:numcoor0,2:5) = tmp_coor;

% compute scaling center of polygonal elements
for iel = 1:numel
    % compute kernel
    kernel = computePolygonKernel( elements{iel}' );
    % compute scaling center
    % compute centroid of the kernel
    [scx,scy] = centroid(polyshape(kernel,'Simplify',false));
    coor(numcoor0+iel,2:5) = [scx,scy,1,1];
end

% coor numbers 
coor(:,1) = (1:numcoor);
% nodes are initially outside
coor(:,7) = -1;
% control point are at the boundary
coor(coor(:,5) == 2, 7) =  0;
 

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
        connectivity{iel}(1,4+n) = nodeNumber;
    end
    % remove repeated values without changing the order
    connectivity{iel}(1,5:end) = unique(connectivity{iel}(1,5:end),'stable');
    % update number of nodes per element
    nel = size(connectivity{iel}(5:end),2); 
    connectivity{iel}(1,4) = nel;
    % maximum number of nodes on any element
    maxnel = max(nel,maxnel);  
    % append correspondig scaling center index
    connectivity{iel}(1,end+1) = numcoor0+iel;  
end

numKnotVectors = countKnotVectors;
idx = 1;
% loop over knotVectors
for ikv = 1:numKnotVectors
    % index initial control point
    knotVectors{ikv}(1,3) = indices(idx);
    % index last control point
    knotVectors{ikv}(1,4) = indices(idx + 1);
    idx = idx + 2;
    
end

numIdxControlPoints = countKnotVectors;
for ictrlp = 1:numIdxControlPoints
    ctrlp = controlPoints_coor{ictrlp};
    nctrlp = idxControlPoints{ictrlp}(1,2);
    
    for n = 1 : nctrlp
    
    a = find ( abs(coor(:,2)-ctrlp(n,1))<1e-10);
    b = find ( abs(coor(:,3)-ctrlp(n,2))<1e-10);
    indices = intersect(a,b);
        
   idxControlPoints{ictrlp}(1,n+2) = indices';
   
    end
    
end

end




