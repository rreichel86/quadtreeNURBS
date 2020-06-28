function [coordinates,element_nodes,nel]=extract_leaf(Quadtree)
%Function extract leaf will help to take out the Quadleaf element data

l = Quadtree.findleaves();
%Gives leaf numbers

[numel]=number_of_elements(Quadtree);
%function to get to know about total number of elements

coordinates = cell(length(l),1);
%This will give an array for storing Quad coordinates

elements = cell(numel,1);
%This will give an array for storing coordinates of elements
element_nodes=cell(numel,1);
j=1;
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
        
        element = [quad,extract_element];
        
        inter_cor=[];
        %No intersectional coordinates
    elseif isempty(Quadtree.Node{l(i),1}{3,1})==1
        %if intersection.horizontal data empty then
        %only intersection.vertical and Quad data taken out
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];
        
        
    elseif isempty(Quadtree.Node{l(i),1}{4,1})==1
        %if intersection.vertical data empty then
        %only intersection.horizontal and Quad data taken out
        
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];
        
        
    else
        %intersection.horizontal,intersection.vertical and Quad data taken
        %out
        inter_cor=[Quadtree.Node{l(i),1}{7,1}(:,1),Quadtree.Node{l(i),1}{7,1}(:,end)];
        %intersectional coordiantes are the 1st and last points of control
        %points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4);
        
        cont_points=Quadtree.Node{l(i),1}{7,1};
        
        coordinate=[quad,inter_cor,cont_points];
        
        element = [quad,extract_element,inter_cor];
        
    end
    
    coordinates{i} = coordinate;
    
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
        elements{j} = element;
        j=j+1;
    else
        %loop for two inter_cor points to know in element definition where they exist
        for k=1:2
            [b] = find ( abs(inter_cor(1,k) - element(1,:))<1e-10);
            [d] = find (abs(inter_cor(2,k) - element(2,:))<1e-10);
            col(:,k)= intersect(b,d);
        end
        if col(1,1) < col(1,2)
            %In Quadtree all control_points are defined in clockwise direction
            %so if col 1 value is less than 2nd these 2 commands work
            elements{j}=[element(:,[1:col(1,1)]),cont_points(:,[2:end-1]),element(:,[col(1,2):end])];
            elements{j+1}=[element(:,[col(1,1):col(1,2)]),fliplr(cont_points(:,[2:end-1]))];
            
        else
            elements{j}=[element(:,[1:col(1,2)]),fliplr(cont_points(:,[2:end-1])),element(:,[col(1,1):end])];
            elements{j+1}=[element(:,[col(1,2):col(1,1)]),cont_points(:,[2:end-1])];
        end
        j=j+2;
    end
end

%Given commands are used to remove the repeated coordinates and arrange
% them by the x coordinte value increasing
coordinates = [coordinates{:}]';
tol=1e-10;
coordinates=uniquetol(coordinates,tol,'ByRows',true);
%A=coordinates';
%qref=reshape(A,[],1);
ncoor = size(coordinates,1) ;
nodes = [1:ncoor]' ;
coordinates = [nodes,coordinates];
figure(2)
plot(coordinates(:,2),coordinates(:,3),'.r') ;
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',8);%to set figure all data in one font
set(gca,'FontSize',10);
text(coordinates(:,2),coordinates(:,3),num2str(nodes));

%Element w.r.t node numbers in cells by using numel
maxnel = 0;
for i=1:numel
    %Following loop for giving nodes to elemnt coordintes by using
    %coordinates row
    
    nel = size(elements{i},2); % number of nodes per element
    for n=1 : nel
        [a] = find ( abs(coordinates(:,2) - elements{i}(1,n))<1e-10);
        [b] = find (abs(coordinates(:,3) - elements{i}(2,n))<1e-10);
        row= intersect(a,b);
        element_nodes{i}(1,n) = row;
    end
    element_nodes{i}= unique(element_nodes{i},'stable');
    %to remove the repeated element_nodes values without changing
    %the order
    maxnel = max(nel,maxnel);  % maximum number of nodes on any element
    element_nodes{i}=[i,nel,element_nodes{i}];
end

end




