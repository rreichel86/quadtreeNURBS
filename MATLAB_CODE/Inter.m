function [Pint,U] = Inter(x1,y1,x2,y2,degree,knots,controlPoints,weights,aux)
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
Pint=[];
U=[];
index=[];
tol = 1e-10;
jk=[];
index=[];
n=length(knots)-degree-2;
num_jk=0;

% First we obtain if there is an intersection between an control poligon
% segment and the quad's edge
for i=1:size(controlPoints,2)-1
    % loop over the control points. Each segment of the control poligon is
    % defined by two consecutive control points.
    Px = controlPoints(1,i);
    Py = controlPoints(2,i);
    if ( i == size(controlPoints,2) )
    else
        Qx = controlPoints(1,i+1);
        Qy = controlPoints(2,i+1);
    end
    %Checks 4 orientations between Pi,Qi,A,B
    OI = orientation([x1 y1],[x2 y2],[Px,Py]);
    OII = orientation([x1 y1],[x2 y2],[Qx,Qy]);
    OIII = orientation([Px,Py],[Qx,Qy],[x1 y1]);
    OIV = orientation([Px,Py],[Qx,Qy],[x2 y2]);
    % If the condition is fulfilled there is a intersection. Obtain the
    % indexes of the control point array that defines the intersecting
    % segement and the number of intersections.
    if OI ~= OII && OIII ~= OIV        
        if OI ~= 0
            index=[index i i+1];
        end       
        if OI == 0
            if i==1
            index=[index i i+1];
            else
                index(end) = i+1;
            end
        end
    end    
end
num_inter=length(index)/2;
% Correcting quad's segment coordinates
if abs(x2-x1) < tol
    y1=y1+aux;y2=y2-aux;
elseif abs(y2-y1) < tol
    x1=x1+aux;x2=x2-aux;
end

% Next step is to obtain the first aproximation of the solution
for i = 1:num_inter
    % We loop over the segments of the control poligon that intersect the
    % quad's edge
    if abs(x2-x1) < tol
        cond3 = x1;
        idx2 = 1;
        
    elseif abs(y2-y1) < tol
        cond3 = y1;
        idx2 = 2;
        
    end
    % findIntrscApprox scans the part of the NURBS controlled by the control
    % points that define the intersected segment. It obtains the first
    % approximation of the solution as output (jj)
    [jj,flag] = findIntrscApprox(cond3,idx2, index((2*i-1):2*i),degree,knots,...
        controlPoints,weights,x1,x2,y1,y2,jk);
    if flag==1
        num_jk = num_jk + 1;
        jk=[jk jj];
    end
end

hold on;
% Based on Steffensen's method obtains the closest solution for a
% given tolerance
for i=1:num_jk
    %loops over the approximated solutions
    a = jk(i);
    hh = 1e-2;
    for ii=1:20
        [f] = curvePoint(n,degree,knots,controlPoints,a,weights);
        O1  = lfunc(f(1),f(2),x1,y1,x2,y2);
        err = abs(lfunc(f(1),f(2),x1,y1,x2,y2));
        if err < tol
            %Avoiding to obtain multiplicities,
            if any(abs(U-a)>tol) || isempty(U)
                if idx2==1 && f(2,1)>y1 && f(2,1)<y2
                    Pint=[Pint f(:)];
                    U=[U a];
                end
                if idx2==2 && f(1,1)>x1 && f(1,1)<x2
                    Pint=[Pint f(:)];
                    U=[U a];
                end
            end
            break
        end
        if (a+hh) > 1
            hh=hh/10  ;
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
                if idx2==1 && f(2,1)>y1 && f(2,1)<y2
                    Pint=[Pint f(:)];
                    U=[U a];
                end
                if idx2==2 && f(1,1)>x1 && f(1,1)<x2
                    Pint=[Pint f(:)];
                    U=[U a];
                end
            end
        end
    end       
end
