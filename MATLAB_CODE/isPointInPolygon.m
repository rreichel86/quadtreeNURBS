function ptInPolygon = isPointInPolygon(polygon, point)
% IsPointInPolygon: Determine if a given point is inside a polygon
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
% ptInPolygon -------------------- gives 1 if point is inside 
%                                        0 if point is at the boundary and 
%                                       -1 otherwise
%
%-------------------------------------------------------------------------

% Compute bounding box that enclosed the polygon
% xmin = min(polygon(:,1));
% ymin = min(polygon(:,2));
xmax = max(polygon(:,1));
ymax = max(polygon(:,2));

% Point in question
q = [point(1) point(2)];
% Arbitrary point outside the polygon 
% For example, use bounding box top right corner 
% and shift it.
r = [xmax+10,ymax+10]; 

% counter number of crossings or intersections
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
    
    % Are q, p_(i) and p_(j) collinear?
    if orientation(q,pi,pj) == 0 % yes
        vi = pi-q;
        vj = pj-q;
        % Is q on the segment p_(i)p_(j)?
        if vi(1) * vj(1) < 0 || vi(2) * vj(2) < 0 || norm(vi) < 1e-10 || norm(vj) < 1e-10 % yes
            % q is at the boundary
            ptInPolygon = 0;
            return
        else % no
            % Does also r lie on the line passing through q, p_(i) and p_(j)?
            if orientation(r,pi,pj) == 0 % yes
                wi = pi-r;
                % Are q and r on the same side?
                %          p_(i)            p_(j)
                %    q & r  *----------------* 
                %           *----------------*  q & r
                % Remark: By definition r can not be on the segment p_(i)p_(j)
                if vi(1) * wi(1) < 0 || vi(2) * wi(2) < 0 % no
                    ph = polygon(a,:); % p_(i-1)
                    pk = polygon(d,:); % p_(j+1)
                    % case (iii) or
                    %
                    %                    q
                    %                    |   
                    %                    .    * ...
                    %                    |   /
                    %                    .  /
                    %                    | /
                    %                    ./
                    %                    * p_(i)
                    %                    |
                    %                    |
                    %                    |
                    %                    |
                    %                    * p_(j)
                    %                    .\
                    %                    | \
                    %                    .  \
                    %                    |   \
                    %                    .    * ...
                    %                    |
                    %                    r
                    %
                    % case (iv)
                    %
                    %                    q
                    %                    |   
                    %                    .    * ...
                    %                    |   /
                    %                    .  /
                    %                    | /
                    %                    ./
                    %                    * p_(i)
                    %                    |
                    %                    |
                    %                    |
                    %                    |
                    %                    * p_(j)
                    %                   /.
                    %                  / |
                    %                 /  .
                    %                /   |
                    %           ... *    .
                    %                    |
                    %                    r
                    %
                    if orientation(ph,q,r)*orientation(pk,q,r) < 0 % case (iv)
                        countIntersection = countIntersection + 1;
                    end
                end
            end
        end
    else % no
        % Are p_(i), q and r collinear?
        if orientation(pi,q,r) == 0 % yes
            vi = q-pi;
            vj = r-pi;
            % Is p_(i) on the segment qr?
            if vi(1) * vj(1) < 0 || vi(2) * vj(2) < 0 || norm(vi) < 1e-10 || norm(vj) < 1e-10 % yes
                ph = polygon(a,:); % p_(i-1)
                % case (i) or
                %
                %                    q
                %                    |   
                %                    .    * ...
                %                    |   /
                %                    .  /
                %                    | /
                %                    ./
                %                    * p_(i)
                %                    .\
                %                    | \
                %                    .  \
                %                    |   \
                %                    .    * ...
                %                    |
                %                    r
                %
                % case (ii)
                %
                %                         * ...
                %                        /
                %                       /
                %                      /
                %                p_(i) /
                %       r - . - . - .*. - . - . - q
                %                     \
                %                      \
                %                       \
                %                        \
                %                         * ...
                %
                if orientation(ph,q,r)*orientation(pj,q,r) < 0 % case (ii)
                    countIntersection = countIntersection + 1;
                end
            end
        else % no
           
            % Do segment qr and p_(i)p_(j) intersect?
            if (orientation(q,pi,pj)*orientation(r,pi,pj) < 0) && ...
                (orientation(pi,q,r)*orientation(pj,q,r) < 0) % yes
                countIntersection = countIntersection + 1;
            end
        end
    end
    
    
end

% if number of intersections even, the point is outside the polygon
% otherwise the point is inside the polygon or on the polygon's boundary
if ( mod(countIntersection,2) == 0 )
   ptInPolygon = -1;
else 
   ptInPolygon = 1;
end
