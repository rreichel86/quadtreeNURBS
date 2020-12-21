function [Q_aux, Quadtree,numInterPoints] = splitting(Quadtree,Q_aux,Quad,...
    NURBS,l,k,i,pos_aux)
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

degree = NURBS.degree;
knots = NURBS.knots;
controlPoints = NURBS.controlPoints;
weights = NURBS.weights;

Q_xmin = Quad(1,1);
Q_xmax = Quad(1,2);
Q_ymin = Quad(2,1);
Q_ymax = Quad(2,3);

% Possible intersections with one of the 4 quad's edges
% Intersection with quad's bottom edge
[Px0, Ux0] = Inter(Q_xmin,Q_ymin,Q_xmax,Q_ymin,degree,knots,controlPoints,weights); 
% Intersection with quad's top edge
[Px1, Ux1] = Inter(Q_xmin,Q_ymax,Q_xmax,Q_ymax,degree,knots,controlPoints,weights); 
% Intersection with quad's left edge
[Py0,Uy0] = Inter(Q_xmin,Q_ymin,Q_xmin,Q_ymax,degree,knots,controlPoints,weights); 
% Intersection with quad's right edge
[Py1,Uy1] = Inter(Q_xmax,Q_ymin,Q_xmax,Q_ymax,degree,knots,controlPoints,weights); 


Px = [Px0, Px1]; % [x1,x2;y1,y2]
Ux = [Ux0, Ux1];
if ~isempty(Px)
    plot(Px(1,:),Px(2,:),'bo', 'LineWidth',1.5);
end
Py = [Py0, Py1]; % [x1,x2;y1,y2]
Uy = [Uy0, Uy1];
if ~isempty(Py)
    plot(Py(1,:),Py(2,:),'bo','LineWidth',1.5);
end

% avoid duplicated knots values
U = [Ux Uy];
U = unique(U);

if any(U == 0)
    if any((U-0.65) > 0)
        U(U == 0)=1;
    end
end

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

refQ = position(l,k,i,pos_aux,Quadtree);
Q_aux = refQ;

% The information of the curve contained in the quad, pointers and auxiliar
% variables are stord at the tree data structure at the node assigned to
% the given quad
[Quadtree] = savetree(Q_aux, Quadtree,k, degree, Px, Py, newKnots, controlPoints, knots_new, weights,Quad);
numInterPoints=length(U);

end


