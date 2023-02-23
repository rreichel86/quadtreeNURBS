function NURBS = CalculateNURBS_2(degree,iKnot,jKnot,knots,controlPoints,weights,numDiv)
% CalculateNURBS: Compute NURBS Curve points 
%
% INPUT:
% degree --------------------- NURBS degree
% knots ---------------------- NURBS knot vector
% controlPoints -------------- NURBS control points 
% weights -------------------- NURBS weights
%
% OUTPUT:
% NURBS Curve points, where 
% 1st column ----------------- x coordinates
% 2nd column ----------------- y coordinates
% 3rd column ----------------- parametric coordinates
%
% -------------------------------------------------------------------------
tol = 1e-10;

if ~exist('numDiv','var')
    numDiv = 100;
end

% number of control points - 1;
n = length(controlPoints)-1;

if abs (iKnot - jKnot) < tol
    subs = iKnot;
else    
    subs = iKnot:(jKnot-iKnot(1))*0.01:jKnot;
end

NURBS = zeros(length(subs),3);

for i = 1:length(subs)
    f = curvePoint(n,degree,knots,controlPoints,subs(i),weights);
    NURBS(i,:) = [f;subs(i)];
end
