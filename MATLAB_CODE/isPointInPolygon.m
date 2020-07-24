function ptInPolygon = isPointInPolygon(polygon, point)
% IsPointInPolygon: Determine if a given point is inside a Polygon
%
% INPUT:
% polygon --------------------- [x_1, y_1;
%                                x_2, y_2;...
%                                x_3, y_3;...
%                                ...
%                                x_n, y_n]
%
% point ----------------------- [ptx, pty]
%
% OUTPUT:
% ptInPolygon -------------------- is 1 if point is inside and 0 otherwise
%
%-------------------------------------------------------------------------

% xmin = polygon(:,1);
% ymin = polygon(:,2);
xmax = max(polygon(:,1));
ymax = max(polygon(:,2));

q = [point(1) point(2)];
r = [xmax+10,ymax+10];

countIntersection = 0;
numvertices = size(polygon,1);

for nv = 1 : numvertices
    
    if nv == numvertices
        a = nv - 1;
        b = nv;
        c = 1;
        d = c + 1;
        
    elseif nv == 1
        a = numvertices;
        b = nv;
        c = nv + 1;
        d = nv + 2;
    else
        a = nv - 1;
        b = nv;
        c = nv + 1;
        d = nv + 2;
    end
    
    pi = polygon(b,:);
    pj = polygon(c,:);
    
    
    if orientation(q,pi,pj) == 0
        vi = pi-q;
        vj = pj-q;
        if vi(1) * vj(1) < 0 || vi(2) * vj(2) < 0 || norm(vi) < 1e-10 || norm(vj) < 1e-10
            countIntersection = countIntersection + 1;
            break;
        else
            if orientation(r,pi,pj) == 0
                wi = pi-r;
                if vi(1) * wi(1) < 0 || vi(2) * wi(2) < 0
                    ph = polygon(a,:);
                    pk = polygon(d,:);
                    if orientation(ph,q,r)*orientation(pk,q,r) < 0
                        countIntersection = countIntersection + 1;
                    end
                end
            end

            
        end
        
    else
        if orientation(pi,q,r) == 0
            vi = q-pi;
            vj = r-pi;
            if vi(1) * vj(1) < 0 || vi(2) * vj(2) < 0 || norm(vi) < 1e-10 || norm(vj) < 1e-10
                ph = polygon(a,:);
                if orientation(ph,q,r)*orientation(pj,q,r) < 0
                    countIntersection = countIntersection + 1;
                end
            end
        else
           
            if (orientation(q,pi,pj)*orientation(r,pi,pj) < 0) && ...
                (orientation(pi,q,r)*orientation(pj,q,r) < 0)
            
                countIntersection = countIntersection + 1;
            end
        end
    end
    
    
end

% if number of intersections even, the point is outside the polygon
% otherwise the point is inside the polygon or on the polygon's boundary
ptInPolygon = mod(countIntersection,2);

end
