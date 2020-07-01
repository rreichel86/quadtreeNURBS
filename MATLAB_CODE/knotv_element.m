function [kv_element,kv_num,maxnk]=knotv_element(Quadtree,connectivity,knot_v,numel,rows);
%function to get knot vector with there connectivity and degree in matrix
%form 
% Input 
% Quadtree data 
% connectivity:connectivity matrix of elements
%knot_v:knot vector of that elements that have NURBS in cell form
%numel:number of elements 
% rows:intersectional coordintes number of that
%leafs that have intersections points
% Output
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

kv_element=cell(numel,1);

kv_num=cell(size(rows,1)/2,1);
j=1;
k=1;
m=1;
maxnk=0;
for i=1:length(l)
    
    %This function will take out the intersection.parametric data from
    %Quadtree
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        
       kv_element{j}=[j,0,connectivity{j,1}(1,2:end)];
       j=j+1;
    else
        kv_element{j}=[j,k,connectivity{j,1}(1,2:end)];
        kv_element{j+1}=[j+1,k,connectivity{j+1,1}(1,2:end)];
        degree=Quadtree.Node{l(i),1}{6,1};
        x=rows([m m+1],1)';
        nk=size(knot_v{i,1}(1,2:end),2);
        kv_num{k}=[k,degree,x,nk,knot_v{i,1}(1,2:end)];
        m=m+2;
        j=j+2;
        k=k+1;
    end
maxnk = max(nk,maxnk);  % maximum number of knot values on any element
end
end