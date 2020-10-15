function [Pint,U] = Inter(x1,y1,x2,y2,degree,knots,controlPoints,weights)
% Inter functions obtains the intersection between a NURBS and a given 
% segment. If it exists, the control poligon will also intersect the
% NURBS. The routine first obtains this segment of the control poligon.
% Afterwards obtains an approximation of the solution: At the end it
% obatins the solution with the Stevenson's method for the given accuracy.
%
% Input:
% x1,y1,x2,y2: definition of the segment
% Definition of the NURBS
% Output:
% Pint: physical coordinates of the intersection point(empty if none)
% U: paramrametrical coordinates of the intersection point(empty if none)

% Initialization of variables
Pint = [];
U = [];
tol = 1e-10;

intrsc_1 = [];
num_intrsc_1 = 0;
intrsc_2 = [];
num_intrsc_2 = 0;
a_array = [];
num_a = 0;

% number of control points / basis functions - 1
n = length(knots)-degree-2;
% ncp number of control points 
ncp = size(controlPoints,2);

% First we obtain if there is an intersection between a control polygon
% and the quad's edge

% a vertical edge
if abs(x2 - x1) < tol
    coorVal = x1;
    coorIdx = 1;
    coorInterval = [y1, y2];
    
% a horizontal edge    
elseif abs(y2 - y1) < tol
    coorVal = y1;
    coorIdx = 2;
    coorInterval = [x1, x2];
end

% Loop over the control points
for i = 1:ncp-1
   
    P = controlPoints(:,i);
    
    if i ~= ncp 
        Q = controlPoints(:,i+1);
    end 
    
    dP = ( P(coorIdx) - coorVal );
    dQ = ( Q(coorIdx) - coorVal );
    
    if dP * dQ < 0 
       intrsc_2 = [intrsc_2 i i+1];
    elseif abs(dP) < tol   
       intrsc_1 =  [intrsc_1 i];
    elseif  abs(dQ) < tol
       intrsc_1 =  [intrsc_1 i+1];
    end     
    
        
end     


if isempty(intrsc_2) && isempty(intrsc_1)
    return
end
if  ~isempty(intrsc_2)
    
    num_intrsc_2 = length(intrsc_2)/2;
end  
if  ~isempty(intrsc_1)  
    
    intrsc_1 = unique(intrsc_1,'stable');
    num_intrsc_1 = length(intrsc_1);
    
end 

% Next step is to obtain an approximation of the solution 
for i = 1:num_intrsc_2
   
    % findIntrscApprox scans the part of the NURBS controlled by the control
    % points that define the intersected segment. It obtains the first
    % approximation of the solution as output (a_pt)
    
    [a_pt,flag] = findIntrscApprox(coorVal,coorIdx,intrsc_2((2*i-1):2*i),degree,knots,...
        controlPoints,weights,x1,x2,y1,y2,a_array);
    if flag == 1
        num_a = num_a + 1;
        a_array = [a_array a_pt];
    end
    
end

for i = 1:num_intrsc_1
   
    % findIntrscApprox scans the part of the NURBS controlled by the control
    % points that define the intersected segment. It obtains the first
    % approximation of the solution as output (a_pt)
    
    [a_pt,flag] = findIntrscApprox(coorVal,coorIdx,intrsc_1(i:i),degree,knots,...
        controlPoints,weights,x1,x2,y1,y2,a_array);
    if flag == 1
        num_a = num_a + 1;
        a_array = [a_array a_pt];
    end
    
end

if (num_a == 0) 
    return
end    

a_knots_array = a_array(3,:);
% avoid duplicated start values
a_knots_array = unique(a_knots_array);
num_a = length(a_knots_array);

% Based on Steffensen's method obtains the closest solution for a
% given tolerance
for i=1:num_a
    %loops over the approximated solutions
    a = a_knots_array(i);
    hh = 1e-2;
    for ii=1:20
        [f] = curvePoint(n,degree,knots,controlPoints,a,weights);
        O1  = lfunc(f(1),f(2),x1,y1,x2,y2);
        err = abs(lfunc(f(1),f(2),x1,y1,x2,y2));
        if err < tol
            %Avoiding to obtain multiplicities,
            if any(abs(U-a)>tol) || isempty(U)
                if coorIdx == 1 
                    if f(2) > y1 && f(2) < y2
                        Pint = [Pint f(:)];
                        U = [U a];
                    elseif ( abs( f(2) - y1 ) < tol ) ||  ( abs( f(2) - y2 ) < tol )
                        Pint = [Pint f(:)];
                        U = [U a];
                    end 
                
                elseif coorIdx == 2
                    if f(1) > x1 && f(1) < x2
                        Pint = [Pint f(:)];
                        U = [U a];
                    elseif ( abs( f(1) - x1 ) < tol ) ||  ( abs( f(1) - x2 ) < tol )
                        Pint = [Pint f(:)];
                        U = [U a];
                    end    
                end
            end
            break
        end    
        if ( abs( a - knots(end) ) < 1e-16 ) ||  ( a > knots(end) )
            break;
        end     
        
        if (a+hh) > knots(end)
            hh=hh/10;
        end
        [g] = curvePoint(n,degree,knots,controlPoints,a+hh,weights);
        [O3] = lfunc(g(1),g(2),x1,y1,x2,y2);
        a =  a - O1/((O3-O1)/(hh));
        if a < 0
            a = abs(a);
        end
        hh = hh/10;
        if ii==20
            %Avoiding to obtain multiplicities
            if any(abs(U-a)>tol) || isempty(U)
                if coorIdx == 1 
                    if f(2) > y1 && f(2) < y2
                        Pint = [Pint f(:)];
                        U = [U a];
                    elseif ( abs( f(2) - y1 ) < tol ) ||  ( abs( f(2) - y2 ) < tol )
                        Pint = [Pint f(:)];
                        U = [U a];
                    end 
                
                elseif coorIdx == 2
                    if f(1) > x1 && f(1) < x2
                        Pint = [Pint f(:)];
                        U = [U a];
                    elseif ( abs( f(1) - x1 ) < tol ) ||  ( abs( f(1) - x2 ) < tol )
                        Pint = [Pint f(:)];
                        U = [U a];
                    end    
                end
            end
        end
    end       
end
