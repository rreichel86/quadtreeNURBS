function [s,flag] = findinterval(cond3,idx, index,degree,knots,...
    controlPoints,weights,x1,x2,y1,y2,jk)                 
% The function obtains the section of the curve where the Control Points 
% that define the intersecting segment have influence, i.e. where their
% associated basis functions are not banishing. Afterwards it obtains the
% first apprimation. This is obtained by subtracting to all the x
% coordinates of the section of the curve, the x corrdinate of a vertical
% segment of the quad, y coordinate for an horizontal quad side. It checks
% where a change of sign of the array occurs, that point of the segment of
% the curve will be used later as the first approximation of the solution
%
% Input: 
% cond3: x/y coordinate of the quad edge
% index: array containing the indixes of the control poligon
% Definition of the NURBS
% Geometrical definition of the quad's segment
% Output:
% s: parametrical coordinate of the first approximation
% flag: flag that indicates if there exists an actual intersection
knotsNew=knots(index(1):(index(end)+degree+1));
n=length(knots)-degree-2;
counter=1;
flag=1;
tol=1e-10;
% If we have a second intersection we start from the other side of the
% curve, that means going backwards through the knot parametrical space
% defined by the knot vector
if isempty(jk)
    for i=knotsNew(1):(knotsNew(end)-knotsNew(1))*0.001:knotsNew(end)
        [f] = curvePoint(n,degree,knots,controlPoints,i,weights);
        P(1:3,counter) = [  f(1);f(2);i];
        counter = counter + 1;
    end
elseif ~isempty(jk)
    for i=knotsNew(end):-(knotsNew(end)-knotsNew(1))*0.001:knotsNew(1)
        [f] = curvePoint(n,degree,knots,controlPoints,i,weights);
        P(1:3,counter) = [  f(1);f(2);i];
        counter = counter + 1;
    end
end



i=1;
g(1)=P(idx,i)-cond3;
while i<length(P)
    g(i+1)=P(idx,i+1)-cond3;
    if idx==1
        %y2>y1 always, from Inter input at reparametrize
        if(g(i)*g(i+1))<=0 && P(2,i)+tol>y1 && P(2,i)-tol<y2 || ...
               ((g(i)*g(i+1))<=0 && P(1,i+1)+tol>y1 && P(1,i+1)-tol<y2)
            break
        end
    end 
    if idx==2 
        %x2>x1 always, from Inter input at reparametrize
        if((g(i)*g(i+1))<=0 && P(1,i)+tol>x1 && P(1,i)-tol<x2) || ...
               ((g(i)*g(i+1))<=0 && P(1,i+1)+tol>x1 && P(1,i+1)-tol<x2)
            break
        end
    end 
    
i=i+1;

end
    %avoiding duplicites
    if i==length(P);flag=0;end
    if ~isempty(jk);if any(abs(jk-P(3,i))<1e-10);flag=0;end;end
    if ~isempty(jk) && i~=1 && any(abs(jk-P(3,i-1))<1e-10);flag=0;end
    if ~isempty(jk) && i~=length(P) && any(abs(jk-P(3,i+1))<1e-10);flag=0;end
    if any(abs(jk)<1e-10) && abs(P(3,i)-1)<1e-10;flag=0;end
s=P(3,i);