function [intrscP,intrscU] = Inter(A,B,intervalTyp,degree,knots,controlPoints,weights)
% Inter: obtain the intersection between a NURBS curve and 
% a given line segment AB.
%
% INPUT:
% A, B ------------- geometrical definition of line segment
% intervalTyp ------ 1 ()
%                    2 (]
%                    3 [)
%                    4 []
% Definition of the NURBS
% degree --------------------- NURBS degree
% knots ---------------------- NURBS knot vector
% controlPoints -------------- NURBS control points 
% weights -------------------- NURBS weights
%
% OUTPUT:
% intrscP: physical coordinates of the intersection point (empty if none)
% intrscU: parametrical coordinates of the intersection point (empty if none)
% 
% -------------------------------------------------------------------------

tol = 1e-12;

% Initialization of variables
intrscP = [];
intrscU = [];
numIntrsc = 0; 

intrsc_1 = [];
num_intrsc_1 = 0;
intrsc_2 = [];
num_intrsc_2 = 0;
StartValues = [];
numStartValues = 0;

% number of control points / basis functions - 1
n = length(knots)-degree-2;
% ncp number of control points 
ncp = size(controlPoints,2);

% First we obtain if there is an intersection between a control polygon
% and the current quad's edge

% Loop over the control points
for i = 1:ncp-1
   
    P = controlPoints(:,i);
    
    if i ~= ncp 
        Q = controlPoints(:,i+1);
    end 
        
    O1 = orientation(P,A,B);
    O2 = orientation(Q,A,B);
    
    if O1*O2 < 0
        intrsc_2 = [intrsc_2 i i+1];
    elseif O1 == 0    
        intrsc_1 =  [intrsc_1 i];
    elseif O2 == 0    
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
    % approximation of the solution as output (u0 )
    
    [u0 ,flag] = findIntrscApprox(intrsc_2((2*i-1):2*i),degree,knots,...
        controlPoints,weights,A,B,StartValues);
    if flag == 1
        numStartValues = numStartValues + 1;
        StartValues = [StartValues u0 ];
    end
end

for i = 1:num_intrsc_1
   
    % findIntrscApprox scans the part of the NURBS controlled by the control
    % points that define the intersected segment. It obtains the first
    % approximation of the solution as output (u0 )
    
    [u0 ,flag] = findIntrscApprox(intrsc_1(i:i),degree,knots,...
        controlPoints,weights,A,B,StartValues);
    if flag == 1
        numStartValues = numStartValues + 1;
        StartValues = [StartValues u0 ];
    end
end

if (numStartValues == 0) 
    return
end    

% Newton-Raphson iteration 

% start values
uStartValues = StartValues(3,:);
% avoid duplicated start values
uStartValues = unique(uStartValues);
% number of start values 
numStartValues = length(uStartValues);


% intrscP = [];
% intrscU = [];

% Line that passes through A and B
tVec = B - A; % tangent vector 
nVec = [-tVec(2);tVec(1)]; % normal vector
nVec = nVec/norm(nVec); % normalized normal vector

for i = 1:numStartValues
    u = uStartValues(i);
    
    for j = 1 : 10
        [R,R_xi] = shape_func_NURBS_1d(u,degree,knots,weights);
        
        C = controlPoints*R;
        dCdu = controlPoints*R_xi;
        
        G = nVec' * ( C - A);
        dG = nVec' * dCdu;
        
        if abs(G) < tol
            
            % Avoiding to obtain multiplicities
            if any(abs(intrscU-u)>tol) || isempty(intrscU)
                 
                alpha = (C - A)'*tVec/(tVec'*tVec);
                
                % check if intersection point lies on line segment AB
                if alpha > 0 && alpha < 1
                    intrscP = [intrscP C];
                    intrscU = [intrscU u];
                elseif abs(alpha) < tol 
                    if intervalTyp == 3 || intervalTyp == 4
                        intrscP = [intrscP C];
                        intrscU = [intrscU u];
                    end 
                elseif abs(alpha - 1) < tol
                   if intervalTyp == 2 || intervalTyp == 4
                        intrscP = [intrscP C];
                        intrscU = [intrscU u];
                    end 
                end                       
            end
            break
        end
        
        u = u - G/dG;

        if ( u > knots(end) )
            u = knots(end);
        elseif ( u < knots(1) )
            u = knots(1); 
        end    
        
    end
    
end

