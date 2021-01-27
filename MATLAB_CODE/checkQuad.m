function nPoints = checkQuad(polygon, points )
% checkQuad: counts the control points in given Quad

counter = 0;
currentPoints = [];
nPoints = size(points,1);


for ii = 1:nPoints
    
      if isPointInPolygon(polygon, points(ii,:))  ~= -1
        
        counter = counter + 1;
        currentPoints(counter,:) = points(ii,:);
    
      end
    
end

% check for dublicated control point 
if ~isempty(currentPoints)
    counter = size(unique(currentPoints,'rows'),1);
end 

nPoints = counter;
end
