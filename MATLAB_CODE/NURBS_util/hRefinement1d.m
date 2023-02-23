
function NURBS = hRefinement1d(NURBS,numRef)
% hRefinement1d: h-refinement for NURBS curves

degree = NURBS.degree;
knots = NURBS.knots;
controlPoints = NURBS.controlPoints;
weights = NURBS.weights;

for i = 1:numRef
    knts = unique(knots);
    insKnots = knts(1:end-1) + 0.5*diff(knts);
    
    for j= 1:length(insKnots)
        nKnt = insKnots(j);
        [knots, controlPoints, weights] = CurveKnotIns(degree,...
            knots, controlPoints, weights, nKnt, 1);
    end
end

NURBS.degree = degree;
NURBS.knots = knots;
NURBS.controlPoints = controlPoints;
NURBS.weights = weights;

end