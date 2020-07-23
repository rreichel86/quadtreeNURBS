function NURBS = CalculateNURBS(degree,knots,controlPoints,weights)
%This function gives as output a NURBS matrix corresponding:
%1st row: x coordinates physical space
%2nd row: y coordinates physical space
%3rd row: parametrical coordinates
%Input: parameters that define the NURBS 

n=length(controlPoints)-1; %number of knot spans
subs = knots(1):(knots(end)-knots(1))*0.001:knots(end);

NURBS = zeros(length(subs),3);

for i = 1:length(subs)
    %loop over parametrical space defined by knot vector
    %obtain each coordinate of the NURBS and store it 
    [f] = curvePoint(n,degree,knots,controlPoints,subs(i),weights);
    
    NURBS(i,:) = [f(1);f(2);subs(i)];
    
    
end