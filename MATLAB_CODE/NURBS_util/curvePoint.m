function C = curvePoint(n,degree,knots,controlPoints,u,weights)
%Output: point C corresponding to the "u" parametrical coordinate
%Input: n,degree,knots,controlPoints,u,weights for the considered NURBS
 
ndim = size(controlPoints,1);
nControlPoints = size(controlPoints,2);
Pw = zeros(ndim + 1, nControlPoints);
%Applying perspective map H, defined at the NURBS Book
for i = 1 : ndim
    Pw(i,:) = controlPoints(i,:) .* weights;
end
Pw(3,:) = weights;

spanIndex = FindSpan(n,degree,u,knots); %Current knot span
N = BasisFun(spanIndex, u, degree, knots);%Value of current basis functions 
Cw = zeros(3,1);

for j = 1:degree+1
        Cw = Cw + (N(j) * Pw(:, spanIndex - degree + j ));
end
Cw = Cw/Cw(3);
C = Cw(1:2);


