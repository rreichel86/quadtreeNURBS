function [Q_aux, Quadtree,numInterPoints] = splitting(Quadtree,Q_aux,Quad,...
    NURBS0,l,k,i,pos_aux)
% splitting: ...
%
% INPUT:
% Quadtree ------------------- Quadtree data structure
% ...
% 
% OUTPUT:
% Quadtree ------------------- Quadtree data structure
% numInterPoints ------------- number of intersection points
% 
% -------------------------------------------------------------------------


Ix = [];
Iy = [];
newKnotVals = [];
numInterPoints = 0;
NURBS_segment = struct([]);

% quad's reference
refQ = position(l,k,i,pos_aux,Quadtree);
Q_aux = refQ;

% get NURBS curve
if k == 1 % first level of decomposition
    data = Quadtree.get(1);
    NURBS = data{3};
else % other levels of decomposition
    refFQ = Q_aux(1:end-2);
    idxFQ = findQ(Quadtree,refFQ);
    data = Quadtree.get(idxFQ);
    NURBS = data{6};
end

if ~isempty(NURBS)
    
    % NURBS curve definition
    degree = NURBS.degree;
    knots = NURBS.knots;
    controlPoints = NURBS.controlPoints;
    weights = NURBS.weights;
    
    % quad's vertices
    Q_xmin = Quad(1,1);
    Q_xmax = Quad(1,2);
    Q_ymin = Quad(2,1);
    Q_ymax = Quad(2,3);
    
    % Possible intersections with one of the 4 quad's edges
    % Intersection with quad's bottom edge
    [Px0, Ux0] = Inter(Q_xmin,Q_ymin,Q_xmax,Q_ymin,degree,knots,controlPoints,weights);
    % Intersection with quad's left edge
    [Py0, Uy0] = Inter(Q_xmin,Q_ymin,Q_xmin,Q_ymax,degree,knots,controlPoints,weights);
    % Intersection with quad's top edge
    [Px1, Ux1] = Inter(Q_xmin,Q_ymax,Q_xmax,Q_ymax,degree,knots,controlPoints,weights);
    % Intersection with quad's right edge
    [Py1, Uy1] = Inter(Q_xmax,Q_ymin,Q_xmax,Q_ymax,degree,knots,controlPoints,weights);
    
    % horizontal intersections
    % physical coordinates
    Px = [Px0, Px1]; % [x1,x2;y1,y2]
    % parametric coordinates
    Ux = [Ux0, Ux1];
    Ix = [Px;Ux];
    % plot horizontal intersections
    if ~isempty(Px)
        plot(Px(1,:),Px(2,:),'bo', 'LineWidth',1.5);
    end
    
    % vertical intersections
    % physical coordinates
    Py = [Py0, Py1]; % [x1,x2;y1,y2]
    % parametric coordinates
    Uy = [Uy0, Uy1];
    Iy = [Py;Uy];
    % plot vertical intersections
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
    
    % Knot insertions:
    % insert knots that correspond to the parametric coordinates of the
    % intersection points
    
    % intersection points in parametric coordinates
    newKnotVals = sort(U);
    % copy knot vector, newKnots is the knot vector after Knot isertion
    newKnots = knots;
    
    if ~isempty(U)
        for j = 1:length(newKnotVals)
            % loop over newKnotVals
            % for knotVal in newknots
            KnotVal = newKnotVals(j);
            numKnotIns = degree-sum(abs(newKnots(:) - KnotVal) < 1e-10);
            % Insert knot numKnotsIns times
            if numKnotIns > 0
                [newKnots, controlPoints, weights] = CurveKnotIns(degree,...
                    newKnots, controlPoints, weights, KnotVal, numKnotIns);
            end
        end
    end
    
    % NURBS curve after knot isertion
    newNURBS.degree = degree;
    newNURBS.controlPoints = controlPoints;
    newNURBS.knots = newKnots;
    newNURBS.weights = weights;
    
    % number of intersections
    numInterPoints = length(newKnotVals);
    
    % Extract a segment of the NURBS curve 
    if numInterPoints == 2
        knotInterval = newKnotVals;
        [NURBS_segment] = extractNURBS_segment(knotInterval, newNURBS);
    elseif  numInterPoints > 2
        knotInterval = [newKnotVals(1), newKnotVals(end)];
        [NURBS_segment] = extractNURBS_segment(knotInterval, newNURBS);
    end
    
end

% store information in the tree data structure in the node assigned to
% the current quad
[Quadtree] = savetree(Q_aux, Quadtree,k, Ix, Iy, newKnotVals, NURBS_segment, Quad);

end


