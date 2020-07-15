function [coor,connectivity,numcoor,maxnel,numel,kv_element,kv_num,maxnk]=extract_leaf(Quadtree)
%Function extract leaf will help to take out the Quadleaf element data
%input
%Quadtree Data
%Output
%coor:coordinates with their weight
%coor=[Number,x-coor,y-coor,weight]
%connectivity:connectivity matrix for elements
%connectivity=[Number,number of coordinates in element,coordinates of element]
%ncoor:number of coordintes
%maxnel:maximum number of nodes on any element
%nmel:number of elements
%kv_element:Vector with knot vector number and connectivity vector that
%leaf 
% kv_element=[Number,Knot_v Number(0 if the quad have no NURBS
%curve),connectivity matrix from 2 column to end] 
% kv_num:Vector having degree of NURBS curve,intersectional coordinates number,
% size of knot vector and knot vector values
%kv_num=[Number,degree of NURBS curve,intersectional coordinates number in clockwise direction,
% size of knot vector and knot vector values]
% maxnk:maximum number of knot values on any element


l = Quadtree.findleaves();
%Gives leaf numbers

[numel]=number_of_elements(Quadtree);
%function to get to know about total number of elements

coor = cell(length(l),1);
%This will give an array for storing Quad coordinates

cp_we = cell(length(l),1);
%Control points and it's associated weights 

knot_v= cell(length(l),1);
%knot vector of the NURBS in quad

elements = cell(numel,1);
%This will give an array for storing coordinates of elements

connectivity=cell(numel,1);
%connectivity matrix of elements

inter_coor= cell(length(l),1);
%Cell for all the intersection coordinates that contain NURBS to in Cell
%form

j=1;
m=1;
for i=1:length(l)
    
    %This function will take out the intersection.parametric data from
    %Quadtree
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    [extract_element] = extracting_element(Quadtree,l,i);
    %Function to get the coordinates of neighboring element that are
    %sharing the boundary
    
    if isempty(intersections) || length(intersections) == 1
        %if intersection data is empty then
        %there will be only Quad vertices no intersection points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        %Quad definition
        
        coordinate = [quad];
        
        cp_w = [];%control points and weight
        
        kv  = [];%knot vector from Quadtree
        
        element = [quad,extract_element];%leaf coordinates from quad 
%         definition and neighboring elemnt quads
        
        inter_cor = [];%No intersectional coordinates
        
    elseif isempty(Quadtree.Node{l(i),1}{3,1})==1
        %if intersection.horizontal data empty then
        %only intersection.vertical and Quad data taken out
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
         %Quad definition
         
        cont_points=Quadtree.Node{l(i),1}{7,1};%control points from Quadtree
        
        weights=Quadtree.Node{l(i),1}{9,1};%weights from Quadtree
        
        cp_w =[cont_points',weights'];%control points and weights
        
        kv  =[m,Quadtree.Node{l(i),1}{8,1}];%knot vector from Quadtree
        
        m=m+1;
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points
        
        
    elseif isempty(Quadtree.Node{l(i),1}{4,1})==1
        %if intersection.vertical data empty then
        %only intersection.horizontal and Quad data taken out
        
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);%Quad definition
        
        cont_points=Quadtree.Node{l(i),1}{7,1};%control points from Quadtree
        
        weights=Quadtree.Node{l(i),1}{9,1};%weights from Quadtree
        
        cp_w =[cont_points',weights'];%control points and weights
        
        kv  =[m,Quadtree.Node{l(i),1}{8,1}];%knot vector from Quadtree
        
        m=m+1;
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points
        
        
    else
        %intersection.horizontal,intersection.vertical and Quad data taken
        %out
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);%Quad definition
        
        cont_points=Quadtree.Node{l(i),1}{7,1};%control points from Quadtree
        
        weights=Quadtree.Node{l(i),1}{9,1};%weights from Quadtree
        
        cp_w =[cont_points',weights'];%control points and weights
        
        kv  =[m,Quadtree.Node{l(i),1}{8,1}];%knot vector from Quadtree
        
        m=m+1;
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];%leaf coordinates from quad 
%         definition,neighboring element quads and intersectional points


         
        
    end
    
    cp_we{i} = cp_w;%control points and weights in one cell
    
    knot_v{i}=kv;%knot vectors in one cell
    
    inter_coor{i}=inter_cor';%intersectional coordiantes in one cell 
%     where x and y coordinates in one row
    
    coor{i} = coordinate;%coordinates in one cell
    
    % Given element commands are to remove the -99 number and arrange the leaf coordinates
    % (including quad,neighboring quad and intersectional coordinates)
    % in the anticlockwise direction starting from left bottom
    element = element(element~=-99);
    ncol = size(element, 1);
    element = reshape(element,[2,ncol/2]);
    x=element(1,:);
    y=element(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);%Gives angle in anticlockwise from [-pi,pi]
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order);
    if a(1,1)>min(a(1,2:end))
        %if there is any point angle which is less than the left bottom
        %point than it should be end point
        element(1,:) =[x_1(1,2:end),x_1(1,1)];
        element(2,:) =[y_1(1,2:end),y_1(1,1)];
    else
        element(1,:) =x_1;
        element(2,:) =y_1;
    end
    %if int_cor empty than there is only one element that was itself leaf
    %otherwise it would be two elements
    if isempty(inter_cor)    
        elements{j} = [element];
        j=j+1;
    else
        %loop for two inter_cor points to know in element definition where they exist
        for k=1:2
            [b] = find ( abs(inter_cor(1,k) - element(1,:))<1e-10);
            [d] = find (abs(inter_cor(2,k) - element(2,:))<1e-10);
            col(:,k)= intersect(b,d);
        end
        %Given loop is only for the control points in Clockwise direction 
        if col(1,1) < col(1,2)
            elements{j}=[element(:,[1:col(1,1)]),cont_points(:,[2:end-1]),element(:,[col(1,2):end])];
            elements{j+1}=[element(:,[col(1,1):col(1,2)]),fliplr(cont_points(:,[2:end-1]))];
            
        else
            elements{j}=[element(:,[1:col(1,2)]),fliplr(cont_points(:,[2:end-1])),element(:,[col(1,1):end])];
            elements{j+1}=[element(:,[col(1,2):col(1,1)]),cont_points(:,[2:end-1])];   
        end
        j=j+2;
    end
end

% control points in matrix form
cp_w = cell2mat(cp_we);
% intersection point coordinates in matrix form
inter_coor = cell2mat(inter_coor);
tol=1e-10;
% remove repeated control points
cp_w = uniquetol(cp_w,tol,'ByRows',true);

% coordinates in matrix form
coor = [coor{:}]';
% remove repeated coordinates and rearrange them by increasing x coordiante
% value
tmp_coor = uniquetol(coor,tol,'ByRows',true);

% number of coordinates without counting scaling centers
numcoor0 = size(tmp_coor,1);
% number of scaling centers = number of elements
numsc = numel;
% update numcoor
numcoor = numcoor0 + numsc;

% prealloc coor matrix 
coor = zeros(numcoor, 4);
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
coor(:,1) = (1:numcoor) ;

% delete tmp_coor
clearvars tmp_coor

%Following loop is to find the control points and assign the weights
for k=1:length(cp_w)
    [a] = find ( abs(coor(:,2)-cp_w(k,1))<1e-10);
    [b] = find ( abs(coor(:,3)-cp_w(k,2))<1e-10);
    row= intersect(a,b);
    if isempty(row)~=1
        coor(row,4)=cp_w(k,3);
    end
end
%Following loop is to find the intersectional coordintes number that will 
% be used to establish knot vector number matrix
rows=zeros(size(inter_coor,1),1);
for k=1:size(inter_coor,1)
    [a] = find ( abs(coor(:,2)-inter_coor(k,1))<1e-10);
    [b] = find ( abs(coor(:,3)-inter_coor(k,2))<1e-10);
    row= intersect(a,b);
   rows(k,1)=row; 
end

figure(2)
plot(coor(:,2),coor(:,3),'.r') ;
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',8);%to set figure all data in one font
set(gca,'FontSize',10);
text(coor(:,2),coor(:,3),num2str(nodes));

%Element w.r.t node numbers in cells by using numel
maxnel = 0;
for iel = 1:numel
    %Following loop for giving nodes to element coordintes by using
    %coordinates row
    
    % number of nodes per element
    nel = size(elements{iel},2); 
    for n = 1 : nel
        [a] = find ( abs(coor(:,2) - elements{iel}(1,n))<1e-10);
        [b] = find (abs(coor(:,3) - elements{iel}(2,n))<1e-10);
        row= intersect(a,b);
        connectivity{iel}(1,n) = row;
    end
    
    % remove repeated values without changing the order
    connectivity{iel} = unique(connectivity{iel},'stable');
    % update number of nodes per element
    nel = size(connectivity{iel},2); 
    
    % maximum number of nodes on any element
    maxnel = max(nel,maxnel);  
    
    connectivity{iel} = [iel,nel,connectivity{iel},numcoor0+iel];
    
end

[kv_element,kv_num,maxnk]=knotv_element(Quadtree,connectivity,knot_v,numel,rows);

end




