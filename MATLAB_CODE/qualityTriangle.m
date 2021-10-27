function alpha = qualityTriangle(A,B,C)
% qualityTriangle: compute alpha-quality coefficient for triangle ABC

AB = B-A;
AC = C-A;
BC = C-B;

alpha = 2*sqrt(3)*det([AB,BC])/((AB'*AB)+(AC'*AC)+(BC'*BC));

end 