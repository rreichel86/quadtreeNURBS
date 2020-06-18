function [coordinates,element_nodes,nel]=extract_leaf(Quadtree);
%Function extract leaf will help to take out the Quadleaf Coordinates data

l = Quadtree.findleaves();
%Gives leaf numbers

[nel]=number_of_elements(Quadtree);

coordinates = cell(length(l),1);
%This will give an array for storing Quad coordinates

elements = cell(nel,1);
%This will give an array for storing coordinates OF elements
element_nodes=cell(nel,1);
j=1;
for i=1:length(l);
    
    %This function will take out the intersection.parametric data from
    %Quadtree
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    [extract_element] = extracting_element(Quadtree,l,i);
        
    
    if isempty(intersections) || length(intersections) == 1;
        %if intersection data is empty then
        %there will be only Quad vertices no intersection points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        coodinate=[quad];
        
        element = [quad,extract_element];
        
        inter_cor=[];
    elseif isempty(Quadtree.Node{l(i),1}{3,1})==1
        %if intersection.horizontal data empty then
        %only intersection.vertical and Quad data taken out
        
        inter_cor=Quadtree.Node{l(i),1}{4,1};
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coodinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];
   
        
    elseif isempty(Quadtree.Node{l(i),1}{4,1})==1
        %if intersection.vertical data empty then
        %only intersection.horizontal and Quad data taken out
        
        inter_cor=Quadtree.Node{l(i),1}{3,1};
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coodinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];

        
    else
        %intersection.horizontal,intersection.vertical and Quad data taken
        %out
        
        inter_cor=[Quadtree.Node{l(i),1}{3,1},Quadtree.Node{l(i),1}{4,1}];
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coodinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];
        
    end
    
    %coordinates data of quad stored in cell array
    
    x=coodinate(1,:);
    y=coodinate(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order);
    if a(1,1)>min(a(1,2:end));
    coodinate(1,:) =[x_1(1,2:end),x_1(1,1)];
    coodinate(2,:) =[y_1(1,2:end),y_1(1,1)];
    else
    coodinate(1,:) =x_1;
    coodinate(2,:) =y_1;
    end
% %   to arrange in anticlockwise from left bottom coordinate
    coordinates{i} = coodinate;
 element = element(find(element~=-99));   
 ncol = size(element, 1);
 element = reshape(element,[2,ncol/2]);
 tol=1e-10;
element=(uniquetol(element',tol,'ByRows',true))';
    x=element(1,:);
    y=element(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order); 
    if a(1,1)>min(a(1,2:end));
    element(1,:) =[x_1(1,2:end),x_1(1,1)];
    element(2,:) =[y_1(1,2:end),y_1(1,1)];
    else
    element(1,:) =x_1;
    element(2,:) =y_1;
    end
    if isempty(inter_cor);
    elements{j} = element;
    j=j+1;
    else
        for k=1:2
        [a,b] = find ( abs(inter_cor(1,k) - element(1,:))<1e-10);
        [c,d] = find (abs(inter_cor(2,k) - element(2,:))<1e-10);
        col(:,k)= intersect(b,d);
        end
       elements{j}=element(:,[1:min(col) max(col):end]);
       elements{j+1}=element(:,min(col):max(col));
       j=j+2;
    end
end
coordinates = [coordinates{:}]';
% %Remove the repeated coordinates
tol=1e-10;
coordinates=uniquetol(coordinates,tol,'ByRows',true);
A=coordinates';
qref=reshape(A,[],1);
nnode = size(coordinates,1) ;
nodes = [1:nnode]' ;
coordinates = [nodes,coordinates];
figure(2)
plot(coordinates(:,2),coordinates(:,3),'.r') ;
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',8);%to set figure all data in one font
set(gca,'FontSize',10);
text(coordinates(:,2),coordinates(:,3),num2str(nodes));

%Element w.r.t node numbers in cells
for i=1:nel;
        for n=1 : size(elements{i},2)
            [a,b] = find ( abs(coordinates(:,2) - elements{i}(1,n))<1e-10);
            [c,d] = find (abs(coordinates(:,3) - elements{i}(2,n))<1e-10);
            row= intersect(a,c);
            element_nodes{i}(1,n) = row;
            end
           element_nodes{i}=[i,size(element_nodes{i},2),element_nodes{i}]; 
end

end




