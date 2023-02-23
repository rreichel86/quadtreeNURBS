
function NURBS = pRefinement1d(NURBS,t)
% pRefinement1d: p-refinement for NURBS curves

if t == 0
    return
end

degree0 = NURBS.degree;
knots0 = NURBS.knots;
controlPoints0 = NURBS.controlPoints;
weights0 = NURBS.weights;

[knots,controlPoints, weights] = DegreeElevateCurve(degree0,knots0,controlPoints0,weights0,t);


NURBS.degree = degree0 + t;
NURBS.knots = knots;
NURBS.controlPoints = controlPoints;
NURBS.weights = weights;

end
