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
flag = 0;
s = -99;
tol = 1e-10;
% If we have a second intersection we start from the other side of the
% curve, that means going backwards through the parametrical space
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

% Line that passes through A and B
tVec = B - A; % tangent vector
nVec = [-tVec(2);tVec(1)]; % normal vector
nVec = nVec/norm(nVec); % normalized normal vector

% compute distances from P to the line that passes through A and B
g = nVec' * ( P(:,1:2)' - A );
pos = 0;
    for i = 1: nkts - 1

        p = P(i,1:2)';
        q = P(i+1,1:2)';

        if g(i)*g(i+1) < 0

            % projection of p onto the line that passes through A and B
            p = p + nVec'*(A - p)* nVec;
            alpha = (p - A)'*tVec/(tVec'*tVec);

            % check if projection lies on line segment AB
            if alpha > 0 && alpha < 1
                pos = i;
                break
            elseif abs(alpha) < tol || abs(alpha - 1) < tol
                pos = i;
                break
            end

            % projection of q onto the line that passes through A and B
            q = q + nVec'*(A - q)* nVec;
            alpha = (q - A)'*tVec/(tVec'*tVec);

            % check if projection lies on line segment AB
            if alpha > 0 && alpha < 1
                pos = i;
                break
            elseif abs(alpha) < tol || abs(alpha - 1) < tol
                pos = i;
                break
            end

        elseif  abs(g(i)) < tol

            % projection of p onto the line that passes through A and B
            p = p + nVec'*(A - p)* nVec;
            alpha = (p - A)'*tVec/(tVec'*tVec);

            % check if projection lies on line segment AB
            if alpha > 0 && alpha < 1
                pos = i;
                break
            elseif abs(alpha) < tol || abs(alpha - 1) < tol
                pos = i;
                break
            end

        elseif abs(g(i+1)) < tol

            % projection of q onto the line that passes through A and B
            q = q + nVec'*(A - q)* nVec;
            alpha = (q - A)'*tVec/(tVec'*tVec);

            % check if projection lies on line segment AB
            if alpha > 0 && alpha < 1
                pos = i;
                break
            elseif abs(alpha) < tol || abs(alpha - 1) < tol
                pos = i;
                break
            end       
        end
    end

if pos ~= 0
    flag = 1;
    s = P(pos,:)';
    
    % avoid duplicated intersection approximations in closed curves
    if ~isempty(intrscArr)
        
        if any( abs( intrscArr(1,:) - P(pos,1) ) < tol  ) && ...
                any( abs( intrscArr(2,:) - P(pos,2) )  < tol)
            
            flag = 0;
            s = -99;
        end
        
    end
    
end


return

