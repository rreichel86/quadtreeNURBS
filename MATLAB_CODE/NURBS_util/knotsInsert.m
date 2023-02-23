function NURBS = knotsInsert(NURBS,insertKnots)
% hRefinement1d: h-refinement for NURBS curves


if isempty(insertKnots)
    return
end


degree = NURBS.degree;
knots = NURBS.knots;
controlPoints = NURBS.controlPoints;
weights = NURBS.weights;

newKnotVals = sort(insertKnots);
% copy knot vector, newKnots is the knot vector after Knot isertion
newKnots = knots;

for j = 1:length(newKnotVals)
    % loop over newKnotVals
    % for knotVal in newknots
    KnotVal = newKnotVals(j);
    numKnotIns = degree-sum(abs(newKnots(:) - KnotVal) < 1e-10);
    % Insert knot numKnotsIns times
    if numKnotIns > 0
        [newKnots, controlPoints, weights] = CurveKnotIns(degree,...
            newKnots, controlPoints, weights, KnotVal, numKnotIns);
    end
end


NURBS.degree = degree;
NURBS.knots = newKnots;
NURBS.controlPoints = controlPoints;
NURBS.weights = weights;

end