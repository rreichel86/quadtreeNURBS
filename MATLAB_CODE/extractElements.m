function [numcoor,coor,numel,connectivity,maxnel,kv_element,kv_num,maxnk]=extractElements(Quadtree)
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
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, nel, node_1,...,node_nel, scaling_center]
% maxnel --------------------- maximum number of nodes on any element
%
%
% -----------------------------------------------------------------------------

% Get Quadtree leaves
leaves = Quadtree.findleaves();
numleaves = length(leaves);

% Get total number of elements
[numel] = countElements(Quadtree);

coor = cell(numleaves,1);
%This will give an array for storing Quad coordinates

controlPoints_coor = cell(numleaves,1);
%Control points and it's associated weights 

knotVectors = cell(numleaves,1);
%knot vector of the NURBS in quad

elements = cell(numel,1);
%This will give an array for storing coordinates of elements

connectivity = cell(numel,1);
%connectivity matrix of elements

intersections_coor = cell(numleaves,1);
%Cell for all the intersection coordinates that contain NURBS to in Cell
%form

j = 1;
countKnotVectors = 0;
% loop over leaves
for i = 1:numleaves
    
    %This function will take out the intersection.parametric data from
    %Quadtree
    intersections = Quadtree.Node{leaves(i),1}{5,1};
    
    % Check neighborhood of current leaf
    % get mid points if they exist
    [midPoints] = extracting_element(Quadtree,leaves,i);

    
    if isempty(intersections) || length(intersections) == 1
        %if intersection data is empty then
        %there will be only Quad vertices no intersection points
        
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
        %Quad definition
        
        coordinates = quad;
        
        cp = [];%control points and weight
        
        kv  = [];%knot vector from Quadtree
        
        element = [quad,midPoints];%leaf coordinates from quad 
%         definition and neighboring elemnt quads
        
        intersectionPoints = [];%No intersection points
        
    elseif isempty(Quadtree.Node{leaves(i),1}{3,1})==1
        %if intersection.horizontal data empty then
        %only intersection.vertical and Quad data taken out
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        %intersection points are the 1st and last points of control
        %points
        
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);
         %Quad definition
         
        controlPoints = Quadtree.Node{leaves(i),1}{7,1}; %control points from Quadtree
        weights = Quadtree.Node{leaves(i),1}{9,1}; %weights from Quadtree
        cp = [controlPoints',weights']; %control points and weights
        
        kv  =[m,Quadtree.Node{leaves(i),1}{8,1}];%knot vector from Quadtree
        countKnotVectors = countKnotVectors + 1;
        
        coordinates = [quad,intersectionPoints,controlPoints];
        
        element = [quad,midPoints,intersectionPoints];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points
        
        
    elseif isempty(Quadtree.Node{leaves(i),1}{4,1})==1
        %if intersection.vertical data empty then
        %only intersection.horizontal and Quad data taken out
        
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        %intersection points are the 1st and last points of control
        %points
        
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);%Quad definition
        
        controlPoints = Quadtree.Node{leaves(i),1}{7,1};%control points from Quadtree
        weights = Quadtree.Node{leaves(i),1}{9,1};%weights from Quadtree
        cp = [controlPoints',weights'];%control points and weights
        
        kv  =[m,Quadtree.Node{leaves(i),1}{8,1}];%knot vector from Quadtree
        countKnotVectors = countKnotVectors + 1;

        coordinates = [quad,intersectionPoints,controlPoints];
        
        element = [quad,midPoints,intersectionPoints];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points
        
        
    else
        %intersection.horizontal,intersection.vertical and Quad data taken
        %out
        intersectionPoints = [Quadtree.Node{leaves(i),1}{7,1}(:,1),Quadtree.Node{leaves(i),1}{7,1}(:,end)];
        %intersection points are the 1st and last points of control
        %points
        
        quad = Quadtree.Node{leaves(i),1}{10,1}(1:2,1:4);%Quad definition
        
        controlPoints = Quadtree.Node{leaves(i),1}{7,1};%control points from Quadtree
        weights = Quadtree.Node{leaves(i),1}{9,1};%weights from Quadtree
        cp = [controlPoints',weights'];%control points and weights
        
        kv  =[m,Quadtree.Node{leaves(i),1}{8,1}];%knot vector from Quadtree
        countKnotVectors = countKnotVectors + 1;
        
        coordinates = [quad,intersectionPoints,controlPoints];
        
        element = [quad,midPoints,intersectionPoints];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points


         
        
    end
    
    controlPoints_coor{i} = cp;%control points and weights in one cell
    
    knotVectors{i}=kv;%knot vectors in one cell
    
    intersections_coor{i} = intersectionPoints';%intersection points in one cell 
%     where x and y coordinates in one row
    
    coor{i} = coordinates;%coordinates in one cell
    
    % Given element commands are to remove the -99 number and arrange the leaf coordinates
    % (including quad,neighboring quad and intersection points)
    % in the anticlockwise direction starting from left bottom
    element = element(element~=-99);
    ncol = size(element, 1);
    element = reshape(element,[2,ncol/2]);
    x = element(1,:);
    y = element(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);%Gives angle in anticlockwise from [-pi,pi]
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order);
    if a(1,1)>min(a(1,2:end))
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
        j = j+1;
    else
        %loop for two intersectionPoints points to know in element definition where they exist
        for k = 1:2
            [b] = find ( abs(intersectionPoints(1,k) - element(1,:))<1e-10);
            [d] = find (abs(intersectionPoints(2,k) - element(2,:))<1e-10);
            col(:,k) = intersect(b,d);
        end
        %Given loop is only for the control points in Clockwise direction 
        if col(1,1) < col(1,2)
            elements{j} = [element(:,[1:col(1,1)]),controlPoints(:,[2:end-1]),element(:,[col(1,2):end])];
            elements{j+1} = [element(:,[col(1,1):col(1,2)]),fliplr(controlPoints(:,[2:end-1]))];
            
            
        else
            elements{j} = [element(:,[1:col(1,2)]),fliplr(controlPoints(:,[2:end-1])),element(:,[col(1,1):end])];
            elements{j+1} = [element(:,[col(1,2):col(1,1)]),controlPoints(:,[2:end-1])];   
        end
        j = j+2;
    end
end

% control points in matrix form
cp = cell2mat(controlPoints_coor);
% intersection point coordinates in matrix form
intersections_coor = cell2mat(intersections_coor);
tol = 1e-10;
% remove repeated control points
cp = uniquetol(cp,tol,'ByRows',true);

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
coor(1:numcoor0,2:3) = tmp_coor;

% compute scaling center of polygonal elements
for iel = 1:numel
    % compute kernel
    kernel = computePolygonKernel( elements{iel}' );
    
    % compute scaling center
    % compute centroid of the kernel
    [scx,scy] = centroid(polyshape(kernel,'Simplify',false));
    coor(numcoor0+iel,2:3) = [scx,scy];
end

% coor numbers 
coor(:,1) = (1:numcoor);
% type node
coor(:,5) = 1;
% node is initially outside
coor(:,7) = -1; 

% delete tmp_coor
clearvars tmp_coor

% loop over control points
for k = 1:length(cp)
    % search control point and assign weight
    [a] = find ( abs(coor(:,2)-cp(k,1))<1e-10);
    [b] = find ( abs(coor(:,3)-cp(k,2))<1e-10);
    row = intersect(a,b);
    if isempty(row)~=1
        coor(row,4) = cp(k,3);
        % type control point
        coor(row,5) = 2;
        % control point is on the boundary
        coor(row,7) = 0;
    end
end


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
        [b] = find (abs(coor(:,3) - elements{iel}(2,n))<1e-10);
        nodeNumber= intersect(a,b);
        connectivity{iel}(1,n) = nodeNumber;
    end
    % remove repeated values without changing the order
    connectivity{iel} = unique(connectivity{iel},'stable');
    % update number of nodes per element
    nel = size(connectivity{iel},2); 
    % maximum number of nodes on any element
    maxnel = max(nel,maxnel);  
    
    connectivity{iel} = [iel,0,nel,connectivity{iel},numcoor0+iel];  
end


end




