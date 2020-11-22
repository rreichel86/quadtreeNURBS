function nPoints = checkQuad( minCoords, maxCoords, Points )
%CHECKQUAD counts the control points in given area
%   The area is defined with the inputs minCoords and maxCoords
%   If a control point (given with the input controlPoints) is between
%   the maximum and minumim coordinates of the area, it can be concluded
%   that this point is in the area, so the counter gets bigger by one for
%   every point in the area. At the end, nPoin gives the total amount of
%   points in the area
counter = 0;
currentPoints = [];
nPoints = size(Points,1);


for ii = 1:nPoints
    
    if  isPointInQuad(minCoords, maxCoords, Points(ii,:)) == 1
        
        counter = counter + 1;
        currentPoints(counter,:) = Points(ii,:);
    
    end
    
end

% Check for dublicated control point 
counter = size(unique(currentPoints,'rows'),1);

nPoints = counter;
end