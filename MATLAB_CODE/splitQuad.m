function [ child1, child2, child3, child4 ] = splitQuad( minX, maxX, minY, maxY )
%SPLITQUAD generates 4 childs form 1 parent
%   The inputs are the min and max coordinates of the parent,
%   the outputs are 4 childs, each with coordinates [minX, maxX ; minY , maxY]

child3 = [(maxX+minX)/2 , maxX ; (maxY+minY)/2 , maxY ];    
child4 = [(maxX+minX)/2 , maxX ; minY , (maxY+minY)/2 ];
child2 = [ minX , (maxX+minX)/2 ; minY , (maxY+minY)/2 ];
child1 = [ minX , (maxX+minX)/2 ; (maxY+minY)/2 , maxY ];

end

