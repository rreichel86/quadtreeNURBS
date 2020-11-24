function s = orientation(Q,R,S)
% orientation: obtain the orientation of three given points Q,R,S 
% s is the output and can take the values: 
%   1 if the given points are oriented CCW 
%  -1 if the given points are oriented CW
%   0 if the given points are collinear

 tol = 1e-10;
 d = (R(1)-Q(1))*(S(2)-Q(2))-(R(2)-Q(2))*(S(1)-Q(1));
 
 if (d < 0) && (abs(d) > tol)
     s = -1;
 elseif d > 0  && (abs(d) > tol)    
     s = 1;
 else
     s = 0;
 end     