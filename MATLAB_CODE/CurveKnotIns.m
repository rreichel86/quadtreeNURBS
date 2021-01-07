function [knots_new, Q, weights] = CurveKnotIns(degree, knots, controlPoints, weights, newKnot, numKnotIns)
% Adapted from the NURBS book (Algorithm A5.1)
p = degree;
np = length(controlPoints);
s = 0;
r = numKnotIns;

k = find(knots >= newKnot(1),1) - 1;

Pw = [controlPoints(1,:) .* weights ;controlPoints(2,:) .* weights; weights];

mp = np + degree + 1;
nq = np + r;

for i = 1:k
    knots_new(i) = knots(i);
end

for i = 1:r
    knots_new(k+i)= newKnot;
end

for i = k+1:mp
    knots_new(i+r) = knots(i);
end

% Save unaltered control points 
for i = 1:k-p
    Qw(:,i) = Pw(:,i);
end

for i = k-s:np
    Qw(:,i+r) = Pw(:,i);
end

for i = 1:p-s+1
    Rw(:,i) = Pw(:,k-p+i-1);
end

% Insert the knot r times
for j = 1:r
    L = k-p+j;
    
    for i = 1:p-j-s+1
        alpha = (newKnot - knots(L+i-1))/(knots(i+k) - knots(L+i-1));
        Rw(:,i) = alpha*Rw(:,i+1) + (1.0-alpha)*Rw(:,i);
    end
    
    Qw(:,L) = Rw(:,1);
    Qw(:,k+r-j-s) = Rw(:,p-j-s+1);
end

% Load remaining control points
for i = (L+1):(k-s-1)
    Qw(:,i) = Rw(:,i-L+1);
end

Q = Qw(1:2,:) ./ Qw(3,:);
weights = Qw(3,:);

end

