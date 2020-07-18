function ptInQuad = isPointInQuad(minCoords,maxCoords,point)
% IsPointInQuad: Determine if a given point is inside a Quad
% 
% INPUT:  
% minCoords ------------------- [xmin, ymin]
% maxCoords ------------------- [xmax, ymax]
% point     ------------------- [ptx, pty]
%
% OUTPUT:
% ptInQuad -------------------- is 1 if point is inside and 0 otherwise
%
%-------------------------------------------------------------------------

xmin = minCoords(1);
ymin = minCoords(2);
xmax = maxCoords(1); 
ymax = maxCoords(2); 

ptx = point(1);
pty = point(2);

ptInQuad = 0;

if (ptx < xmin) 
    return
elseif (ptx > xmax) 
    return
elseif (pty < ymin) 
    return
elseif (pty > ymax) 
    return
else
    ptInQuad = 1;
end

end