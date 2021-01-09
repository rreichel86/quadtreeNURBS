function [Pint,U] = Inter(x1,y1,x2,y2,degree,knots,controlPoints,weights)
% Inter: obtain the intersection between a NURBS curve and a given line segment. 
%
% INPUT:
% x1, y1, x2, y2 ------------- geometrical definition of line segment
% Definition of the NURBS
% degree --------------------- NURBS degree
% knots ---------------------- NURBS knot vector
% controlPoints -------------- NURBS control points 
% weights -------------------- NURBS weights
%
% OUTPUT:
% Pint: physical coordinates of the intersection point (empty if none)
% U: parametrical coordinates of the intersection point (empty if none)
% 
% -------------------------------------------------------------------------

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
% number of start values 
num_a = length(a_knots_array);


p = [x1;y1];
q = [x2;y2];
tVec = q - p;
nVec = [-tVec(2);tVec(1)];
nVec = nVec/norm(nVec);


for i = 1:num_a
    a = a_knots_array(i);
    
    for j = 1 : 10
        [R,R_xi] = shape_func_NURBS_1d(a,degree,knots,weights);
        
        f = controlPoints*R;
        df = controlPoints*R_xi;
        
        G = nVec' * ( f - p);
        dG = nVec' * df;
        
        if abs(G) < tol
%         if abs(f(coorIdx) - coorVal) < tol
            
            % Avoiding to obtain multiplicities
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
        
        a = a - G/dG;
%         a = a - (f(coorIdx) - coorVal)/df(coorIdx);
        
        
        if ( a > knots(end) )
            a = knots(end);
        elseif ( a < knots(1) )
            a = knots(1); 
        end    
        
    end
    
end

