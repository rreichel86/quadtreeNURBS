function OP = NURBSCurveOrientation(NURBS)

ncp = length(NURBS.knots) - 1 -NURBS.degree;
poly = NURBS.controlPoints.';

% determine signed polygon area
% area =  1 for CCW
% area = -1 for CW
area = poly(ncp,1) * poly(1,2) - poly(1,1) * poly(ncp,2);

for i = 1: ncp-1
    area = area + poly(i,1) * poly(i+1,2);
    area = area - poly(i+1,1) * poly(i,2);
end

% - area is positive, thus vertices are arranged CCW, OP = 1
% - area is negative, thus vertices are arranged CW, OP = -1 

if area > 0
  OP = 1;
else
  OP = -1;
end


end