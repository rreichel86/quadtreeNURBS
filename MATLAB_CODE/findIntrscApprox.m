function [s,flag] = findIntrscApprox(coorVal,coorIdx,interval,degree,knots,...
    controlPoints,weights,x1,x2,y1,y2,intrscArr)
% findIntrscApprox: determine first intersection approximation between
% a Quad edge and a NURBS curve.
%
% INPUT:
% coorVal ----------------------- x or y coordinate of Quad edge
% coorIdx ------------------------ index of the corresponding Quad edge coordiante
% interval ------------------- knot vector interval 
% degree --------------------- NURBS degree
% knots ---------------------- NURBS knot vector
% controlPoints -------------- NURBS control points 
% weights -------------------- NURBS weights
% x1, y1, x2, y2 ------------- geometrical definition of Quad edge
% intrscArr ------------------ record of previous intersections
%
% OUTPUT:
% s -------------------------- first intersection approximation 
%                              in parametric coordinates
% flag ----------------------- flag that indicates if there exists 
%                              an actual intersection
%
% -------------------------------------------------------------------------

kN = knots(interval(1):(interval(end)+degree+1));
flag = 1;
tol = 1e-10;
% If we have a second intersection we start from the other side of the
% curve, that means going backwards through the knot parametrical space
if isempty(intrscArr)
    kts = kN(1):(kN(end)-kN(1))*0.01:kN(end);
    nkts = length(kts);
elseif ~isempty(intrscArr)
    kts = kN(end):-(kN(end)-kN(1))*0.01:kN(1);
    nkts = length(kts);
end

% P = [x, y, s]
P = zeros(nkts,3);
n = length(knots)-degree-2;
counter = 0;
for k = kts
    point = curvePoint(n,degree,knots,controlPoints,k,weights);
    counter = counter + 1;
    P(counter,:) = [point' k];
    
end

        end
    end
        end
end
    end
    
end

