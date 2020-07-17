function [Q_aux, Quadtree,numInterPoints] = splitting(Quadtree, Q_aux, Quad, ...
        controlPoints, knots, weights, degree,l, k, i, pos_aux)
% Splitting function obtains the description of an eventual segment of the 
% NURBS contained in the quad.First we need to detect the intersection 
% points, then performorm knot insertion for the intersection untill C0 
% continuity. At the mathematical description of the segment of the NURBS, 
% pointers and other variables are stored the tree data structure.  

% Input: 
% Quadtree data
% Quad: geometrical definition of the considered quad 
% Definition of the NURBS
% Auxiliar variables previously defined
% Output:
% Quadtree data structure cointining the splitted segment of the NURBS


% First we obtain the possible intersections with the four quad's edges
% Intersection down horizontal
aux=0.9*max(abs(Quad(1,:)));
[Px0, Ux0] = Inter(min(Quad(1,:))-aux,min(Quad(2,:)),max(Quad(1,:))+aux,...
    min(Quad(2,:)),degree,knots,controlPoints,weights,aux); 
%Intersection up horizontal
[Px1, Ux1] = Inter(min(Quad(1,:))-aux,max(Quad(2,:)),max(Quad(1,:))+aux,...
    max(Quad(2,:)),degree,knots,controlPoints,weights,aux); 
%Intersection left vertical
[Py0,Uy0] = Inter(min(Quad(1,:)),min(Quad(2,:))-aux,min(Quad(1,:)),...
    max(Quad(2,:))+aux,degree,knots,controlPoints,weights,aux); 
%Intersection right vertical
[Py1,Uy1] = Inter(max(Quad(1,:)),min(Quad(2,:))-aux,max(Quad(1,:))...
    ,max(Quad(2,:))+aux,degree,knots,controlPoints,weights,aux); 
Px = [Px0, Px1]; %[x1,x2;y1,y2]
Ux = [Ux0, Ux1];
if ~isempty(Px);plot(Px(1,:),Px(2,:),'bo', 'LineWidth',1.5);end
Py = [Py0, Py1]; %[x1,x2;y1,y2]
Uy = [Uy0, Uy1];
if ~isempty(Py);plot(Py(1,:),Py(2,:),'bo','LineWidth',1.5);end

U = [Ux Uy];
if any(U==0);if any((U-0.65)>0);U(U==0)=1;end;end

% After determining the intersection points we perform a knot insertion,
% beeing the knots the parametrical coordinates of the intersection points
newKnots=sort(U); 
knots_new=knots;
if ~isempty(U)
    for j= 1:length(newKnots)
        % loop over newKnots vector. We insert them at knot vector
        newKnot=newKnots(j);
        num_ins = degree-sum(abs(knots_new(:) - newKnot) < 1e-10);
        if num_ins > 0 
        for ij=1:num_ins % Insert knot untill C0 continuity condition
            [knots_new, controlPoints, weights] = CurveKnotIns(degree,...
                controlPoints, knots_new, weights, newKnot);
        end
        end
    end
end

position_Q = position(l,k,i,pos_aux,Quadtree);
Q_aux = position_Q;

% The information of the curve contained in the quad, pointers and auxiliar
% variables are stord at the tree data structure at the node assigned to
% the given quad
[Quadtree] = savetree(Q_aux, Quadtree,k, degree, Px, Py, newKnots, controlPoints, knots_new, weights,Quad);
numInterPoints=length(U);

end


