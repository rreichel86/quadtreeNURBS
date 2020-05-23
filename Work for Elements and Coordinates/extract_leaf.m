function [cordinates]=extract_leaf(Quadtree);
%Function extract leaf will help to take out the Quadleaf Coordinates data

l = Quadtree.findleaves();
%Gives leaf numbers

cor_1 = cell(length(l),1);
%This will give an array for storing Quad coordinates

for i=1:length(l);
    
    %This function will take out the intersection.parametric data from
    %Quadtree
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        %if intersection data is empty then
        %there will be only Quad vertices no intersection points
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4)
        
        cor=[quad]
        
        
    elseif isempty(Quadtree.Node{l(i),1}{3,1})==1
        %if intersection.horizontal data empty then
        %only intersection.vertical and Quad data taken out
        
        int_ver=Quadtree.Node{l(i),1}{4,1}
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4)
        
        cor=[quad,int_ver]
        
    elseif isempty(Quadtree.Node{l(i),1}{4,1})==1
        %if intersection.vertical data empty then
        %only intersection.horizontal and Quad data taken out
        
        int_hor=Quadtree.Node{l(i),1}{3,1}
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4)
        
        cor=[quad,int_hor]
        
        
    else
        %intersection.horizontal,intersection.vertical and Quad data taken
        %out
        
        int_hor=Quadtree.Node{l(i),1}{3,1}
        
        int_ver=Quadtree.Node{l(i),1}{4,1}
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4)
        
        cor=[quad,int_hor,int_ver]
        
    end
    
    cor_1{i} = cor;
    %coordinates data of quad stored in cell array
    
    x=cor_1{i,1}(1,:);
    y=cor_1{i,1}(2,:);
    cx = mean(x);
    cy = mean(y);
    a = atan2(y - cy, x - cx);
    [~, order] = sort(a, 'ascend');
    x_1 = x(order);
    y_1 = y(order);
    cor_1{i,1}(1,:) =x_1;
    cor_1{i,1}(2,:) =y_1;
    %to arrange in anticlockwise from left bottom coordinate
end
cor = [cor_1{:}]'

[UniXY,Index]=unique(cor,'rows','first')
DupIndex=setdiff(1:size(cor,1),Index)
cor(DupIndex,:)=[];
%Remove the repeated coordinates
cordinates=cor
end




