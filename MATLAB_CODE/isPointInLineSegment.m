function [pointInSegment, locPoint] = isPointInLineSegment(A,B,P,intervalTyp)
% IsPointInLineSegment: Determine if given point P lies on
% line segment AB
%
% INPUT:
% A, B ------------- geometrical definition of line segment
% P ---------------- given point 
% intervalTyp ------ 1 ()
%                    2 (]
%                    3 [)
%                    4 []
%
% OUTPUT:
% pointInSection ---- is 1 if point lies on line segment AB
%                     and 0 otherwise
% locPoint ---------- take a value between 0 and 1 if P lies on line segment AB
%                     and is -99 otherweise
%
%-------------------------------------------------------------------------

tol = 1e-12;

pointInSegment = 0;
locPoint = -99;

% Line that passes through A and B
tVec = B - A; % tangent vector
nVec = [-tVec(2);tVec(1)]; % normal vector
nVec = nVec/norm(nVec); % normalized normal vector

% check if P lies on line that passes through A and B
if nVec' * (P - A)  < tol
    
    alpha = (P - A)'*tVec/(tVec'*tVec);
    
    % check if P lies on line segment AB
    if alpha > 0 && alpha < 1
        pointInSegment = 1;
        locPoint = alpha;
    elseif abs(alpha) < tol
        if intervalTyp == 3 || intervalTyp == 4
            pointInSegment = 1;
            locPoint = alpha;
        end
    elseif abs(alpha - 1) < tol
        if intervalTyp == 2 || intervalTyp == 4
            pointInSegment = 1;
            locPoint = alpha;
        end
    end
     
end

end