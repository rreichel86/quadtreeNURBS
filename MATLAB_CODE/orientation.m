function s = orientation(Q,R,S)
% orientation function obtains the orientation of three given points Q,R,S 
% by terms of the cross product. s is the outut and can take the values 1
% for a positive orientation (right-hand rule), -1 for negative, 0 for
% colinear
 tol = 1e-10;
 d = (R(1)-Q(1))*(S(2)-Q(2))-(R(2)-Q(2))*(S(1)-Q(1));
 
 if (d < 0) && (abs(d) > tol)
     s = -1;
 elseif d > 0  && (abs(d) > tol)    
     s = 1;
 else
     s = 0;
 end     