function [s,flag] = findIntrscApprox(coor,idx,interval,degree,knots,...
    controlPoints,weights,x1,x2,y1,y2,intrscArr)
% findIntrscApprox: determine first intersection approximation between
% a Quad edge and a NURBS curve.
%
% INPUT:
% coor ----------------------- x or y coordinate of Quad edge
% idx ------------------------ index of the corresponding Quad edge coordiante
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
    kts = kN(1):(kN(end)-kN(1))*0.001:kN(end);
    nkts = length(kts);
elseif ~isempty(intrscArr)
    kts = kN(end):-(kN(end)-kN(1))*0.001:kN(1);
    nkts = length(kts);
end

% P = [x, y, s]
P = zeros(nkts,3);
n = length(knots)-degree-2;
counter = 1;
for k = kts
    point = curvePoint(n,degree,knots,controlPoints,k,weights);
    P(counter,:) = [point' k];
    counter = counter + 1;
end

g = P(:,idx)-coor;
for i = 1: (length(P)-1)
    if idx == 1
        % y2 > y1 always, from Inter input at reparametrize
        if(g(i)*g(i+1))<=0 && P(i,2)+tol>y1 && P(i,2)-tol<y2 || ...
                ((g(i)*g(i+1))<=0 && P(i+1,2)+tol>y1 && P(i+1,2)-tol<y2)
            break
        end
    end
    if idx == 2
        % x2 > x1 always, from Inter input at reparametrize
        if((g(i)*g(i+1))<=0 && P(i,1)+tol>x1 && P(i,1)-tol<x2) || ...
                ((g(i)*g(i+1))<=0 && P(i+1,1)+tol>x1 && P(i+1,1)-tol<x2)
            break
        end
    end
end

% avoiding duplicites
if i == (length(P) - 1)
    flag = 0;
end
if ~isempty(intrscArr)
    if any(abs(intrscArr-P(i,3))<1e-10)
        flag = 0;
    end
end
if ~isempty(intrscArr) && i~=1 && any(abs(intrscArr-P(i-1,3))<1e-10)
    flag = 0;
end
if ~isempty(intrscArr) && i~=length(P) && any(abs(intrscArr-P(i+1,3))<1e-10)
    flag = 0;
end
if any(abs(intrscArr)<1e-10) && abs(P(i,3)-1)<1e-10
    flag = 0;
end

s = P(i,3);
