function NURBS_pts = CalculateNURBS(NURBS)
% CalculateNURBS: Compute NURBS Curve points
%
% INPUT:
% NURBS.degree --------------------- NURBS degree
% NURBS.knots ---------------------- NURBS knot vector
% NURBS.controlPoints -------------- NURBS control points
% NURBS.weights -------------------- NURBS weights
%
% OUTPUT:
% NURBS Curve points, where
% 1st column ----------------- x coordinates
% 2nd column ----------------- y coordinates
% 3rd column ----------------- parametric coordinates
%
% -------------------------------------------------------------------------

degree = NURBS.degree;
knots = NURBS.knots;
controlPoints = NURBS.controlPoints;
weights = NURBS.weigths;

n = length(controlPoints)-1;
subs = knots(1):(knots(end)-knots(1))*0.01:knots(end);
NURBS_pts = zeros(length(subs),3);

for i = 1:length(subs)
    f = curvePoint(n,degree,knots,controlPoints,subs(i),weights);
    NURBS_pts(i,:) = [f;subs(i)];
end
