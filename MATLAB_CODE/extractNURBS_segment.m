function [NURBS_segment] = extractNURBS_segment(knotVals, NURBS)

% NURBS segment
NURBS_segment = struct;

degree = NURBS.degree;
knots = NURBS.knots;
controlPoints = NURBS.controlPoints;
weights = NURBS.weights;

% extact NURBS segment
k0 = find( abs(knots - knotVals(1) ) < 1e-10, 1, 'First');
k1 = find( abs(knots - knotVals(2) ) < 1e-10, 1, 'Last');

if k0 == 1 && k1 ~= length(knots)
    NURBS_segment.degree = degree;
    NURBS_segment.controlPoints = controlPoints(1:2,k0:(k1-degree));
    NURBS_segment.knots = [knots(k0:k1) knots(k1)];
    NURBS_segment.weights = weights(k0:(k1-degree));
elseif k0 ~= 1 && k1 == length(knots)
    NURBS_segment.degree = degree;
    NURBS_segment.controlPoints = controlPoints(1:2,(k0-1):(k1-degree-1));
    NURBS_segment.knots = [knots(k0) knots(k0:k1)];
    NURBS_segment.weights = weights((k0-1):(k1-degree-1));
elseif k0 ~= 1 && k1 ~= length(knots)
    NURBS_segment.degree = degree;
    NURBS_segment.controlPoints = controlPoints(1:2,(k0-1):(k1-degree));
    NURBS_segment.knots = [knots(k0) knots(k0:k1) knots(k1)];
    NURBS_segment.weights = weights((k0-1):(k1-degree));
else
    NURBS_segment.degree = degree;
    NURBS_segment.controlPoints = controlPoints;
    NURBS_segment.knots = knots;
    NURBS_segment.weights = weights;
end

end