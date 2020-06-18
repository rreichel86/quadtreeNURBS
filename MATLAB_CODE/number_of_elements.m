function [nel]=number_of_elements(Quadtree);

%Count number of elements

l = Quadtree.findleaves();
nel=0;

for i=1:length(l);
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        nel=nel+1;
    else
        nel=nel+2;
    end
end
end