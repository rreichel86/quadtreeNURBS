function [ child1, child2, child3, child4 ] = splitQuad( minX, maxX, minY, maxY )
% splitQuad: generate 4 children from the current parent Quad 
% The inputs are the min and max coordinates of the current parent Quad.
% The outputs are its 4 children, each with coordinates [minX, maxX ; minY , maxY]

child3 = [(maxX+minX)/2 , maxX ; (maxY+minY)/2 , maxY ];    
child4 = [(maxX+minX)/2 , maxX ; minY , (maxY+minY)/2 ];
child2 = [ minX , (maxX+minX)/2 ; minY , (maxY+minY)/2 ];
child1 = [ minX , (maxX+minX)/2 ; (maxY+minY)/2 , maxY ];

end

