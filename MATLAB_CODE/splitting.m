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

tol = 1e-10;
numIntersections = 0;
% array that contains intersection points
Ip = []; 
% array that contains corresponding intersection points location
Is = [];
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
    intersections = data{5};
    numIntersections = length(intersections);
    NURBS = data{6};
end

if  isempty(NURBS) && numIntersections == 1
    
    vertexNum = [4,1,3,2];
    intersectionPoint = data{3};
   
    if norm(Quad(:,vertexNum(i)) - intersectionPoint(1:2)) < tol
      Ip = data{3};
      newKnotVals = data{5}; 
    end     
    
end     

if ~isempty(NURBS)
    
    % NURBS curve definition
    degree = NURBS.degree;
    knots = NURBS.knots;
    controlPoints = NURBS.controlPoints;
    weights = NURBS.weights;
    
    % quad's vertices
    V1 = Quad(:,1);
    V2 = Quad(:,2);
    V3 = Quad(:,3);
    V4 = Quad(:,4);
    
    % Possible intersections with one of the 4 quad's edges

    % Intersection with quad's bottom edge [V1 V2)
    [Px0, Ux0, nIx0] = Inter(V1,V2,3,degree,knots,controlPoints,weights);
    % Intersection with quad's right edge [V2 V3)
    [Py0, Uy0, nIy0] = Inter(V2,V3,3,degree,knots,controlPoints,weights);
    % Intersection with quad's top edge [V3 V4)
    [Px1, Ux1, nIx1] = Inter(V3,V4,3,degree,knots,controlPoints,weights);
    % Intersection with quad's left edge [V4 V1)
    [Py1, Uy1, nIy1] = Inter(V4,V1,3,degree,knots,controlPoints,weights);
    
    % intersection points (physical space)
    Ip = [Px0, Py0, Px1, Py1;...
          Ux0, Uy0, Ux1, Uy1];
      
    % and the corresponding knots
    U = [Ux0, Uy0, Ux1, Uy1];
    % avoid duplicated knot values
    U = unique(U);
    
    if any(U == 0)
        if any((U-0.65) > 0)
            U(U == 0)=1;
        end
    end
    
    % plot intersection points 
    if ~isempty(Ip)
        plot(Ip(1,:),Ip(2,:),'bo', 'LineWidth',1.5);
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
        
         if nIx0 == 1
             Is = [Is 1];
         elseif  nIx0 == 2
             numInterPoints = numInterPoints + 1;
         end  
         
         if nIy0 == 1
             Is = [Is 2]; 
         elseif  nIy0 == 2
             numInterPoints = numInterPoints + 1;
         end
         
         if nIx1 == 1 
             Is = [Is 3];
         elseif  nIx1 == 2
             numInterPoints = numInterPoints + 1;
         end 
         
         if nIy1 == 1 
             Is = [Is 4];
         elseif  nIy1 == 2
             numInterPoints = numInterPoints + 1;
         end 
        
    elseif  numInterPoints > 2
        knotInterval = [newKnotVals(1), newKnotVals(end)];
        [NURBS_segment] = extractNURBS_segment(knotInterval, newNURBS);
    end
    
end

% store information in the tree data structure in the node assigned to
% the current quad
[Quadtree] = savetree(Q_aux, Quadtree,k, Ip, Is, newKnotVals, NURBS_segment, Quad);

end


