function NURBS = CalculateNURBS(degree,knots,controlPoints,weights)
%This function gives as output a NURBS matrix corresponding:
%1st row: x coordinates physical space
%2nd row: y coordinates physical space
%3rd row: parametrical coordinates
%Input: parameters that define the NURBS 

counter = 1;%initialize counter
n=length(controlPoints)-1; %number of knot spans
subs = knots(1):(knots(end)-knots(1))*0.01:knots(end);

NURBS = zeros(3,size(subs,2));

for i = subs
    %loop over parametrical space defined by knot vector
    %obtain each coordinate of the NURBS and store it 
    [f] = curvePoint(n,degree,knots,controlPoints,i,weights);
    
    NURBS(1:3,counter) = [f(1);f(2);i];
    counter=counter+1;
    
end