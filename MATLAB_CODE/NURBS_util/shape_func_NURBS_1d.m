
function [R,R_xi] = shape_func_NURBS_1d(u,degree,knots,weights)

n = size(knots,2) - 1 - degree - 1;

R = zeros(n+1,1);
R_xi = zeros(n+1,1);

spanIndex = FindSpan(n,degree,u,knots);
N = Der1BasisFun(spanIndex,u,degree,knots)';

W = 0;
W_xi = 0;
for i = 1:degree+1
    wgt = weights(spanIndex - degree + i);
    W = W + N(i,1)*wgt;
    W_xi = W_xi + N(i,2)*wgt;
end 

for i = 1:degree+1
    R(spanIndex - degree + i) = N(i,1)*weights(spanIndex - degree + i)/W;
    R_xi(spanIndex - degree + i) = (N(i,2)*W - N(i,1)*W_xi)*weights(spanIndex - degree + i)/(W*W);
end 
