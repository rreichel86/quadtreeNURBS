function [nel]=number_of_elements(Quadtree)
%function to get to know about total number of elements

l = Quadtree.findleaves();
%Calling each leaf of Quadtree
nel=0;

%Following loop is basically if ther are intersectional points than there
%are 2 elements otherwise 1
for i=1:length(l)
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        nel=nel+1;
    else
        nel=nel+2;
    end
end
end