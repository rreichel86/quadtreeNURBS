function nPoints = checkQuad( minCoords, maxCoords, minCoordsGlobal, controlPoints )
%CHECKQUAD counts the control points in given area
%   The area is defined with the inputs minCoords and maxCoords
%   If a control point (given with the input controlPoints) is between
%   the maximum and minumim coordinates of the area, it can be concluded
%   that this point is in the area, so the counter gets bigger by one for
%   every point in the area. At the end, nPoin gives the total amount of
%   points in the area
counter = 0;
currentPoint = [];
for ii = 1:size(controlPoints,2)
    if minCoords(1) == minCoordsGlobal(1) && minCoords(2) == minCoordsGlobal(2)
        if (controlPoints(1,ii) >= minCoords(1) && ...
                controlPoints(1,ii) <= maxCoords(1) ) && ...
                (controlPoints(2,ii) >= minCoords(2) && ...
                controlPoints(2,ii) <= maxCoords(2))
            
            % Makes sure it does not count duplicate control points
            currentPoint(1:2,(size(currentPoint,2)+1)) = controlPoints(:,ii);
            
            counter = size(unique(currentPoint','rows'),1);
            
        end
        nPoints = counter;
    elseif minCoords(1) == minCoordsGlobal(1) && minCoords(2) ~= minCoordsGlobal(2)
        if (controlPoints(1,ii) >= minCoords(1) && ...
                controlPoints(1,ii) <= maxCoords(1) ) && ...
                (controlPoints(2,ii) > minCoords(2) && ...
                controlPoints(2,ii) <= maxCoords(2))
            
            % Makes sure it does not count duplicate control points
            currentPoint(1:2,(size(currentPoint,2)+1)) = controlPoints(:,ii);
            
            counter = size(unique(currentPoint','rows'),1);
        end
    elseif minCoords(1) ~= minCoordsGlobal(1) && minCoords(2) == minCoordsGlobal(2)
        if (controlPoints(1,ii) > minCoords(1) && ...
                controlPoints(1,ii) <= maxCoords(1) ) && ...
                (controlPoints(2,ii) >= minCoords(2) && ...
                controlPoints(2,ii) <= maxCoords(2))
            
            % Makes sure it does not count duplicate control points
            currentPoint(1:2,(size(currentPoint,2)+1)) = controlPoints(:,ii);
            
            counter = size(unique(currentPoint','rows'),1);
        end
    elseif minCoords(1) ~= minCoordsGlobal(1) && minCoords(2) ~= minCoordsGlobal(2)
        if (controlPoints(1,ii) > minCoords(1) && ...
                controlPoints(1,ii) <= maxCoords(1) ) && ...
                (controlPoints(2,ii) > minCoords(2) && ...
                controlPoints(2,ii) <= maxCoords(2))
            
            % Makes sure it does not count duplicate control points
            currentPoint(1:2,(size(currentPoint,2)+1)) = controlPoints(:,ii);
            
            counter = size(unique(currentPoint','rows'),1);
            
        end
        
    end
    
end
nPoints = counter;
end