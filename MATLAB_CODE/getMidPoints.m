function [numMidPoints,midPoints,locMidPoints] = getMidPoints(Quadtree,idxQ,refQ)
% getMidPoints: get current Quad mid points

% current Quad level
levelQ = length(refQ);
% mid points coordinates
midPoints = zeros(2,4) - 99;
% mid points location 
locMidPoints = zeros(1,4);
% number of mid points 
numMidPoints = 0;

% search for current Quad neighbours
% Loop over directions:
% 1 - West
% 2 - South
% 3 - East
% 4 - North
for dir = 1:4
    % determine possible neighbour Quad in current direction
    [exist_NQ, refNQ] = refNeighbour(refQ,dir);
    if exist_NQ == 1
        % look for neighbour Quad reference in reference array
        % and get its position or its ancestor position
        
        idxNQ = findNeighbour(Quadtree,idxQ,refQ, refNQ);
        levelNQ = length( Quadtree.Node{idxNQ,1}{2,1}(1:end) );
        
        % are both of the same size ?
        if levelQ == levelNQ
            
            % has neighbour Quad children ?
            %children = Quadtree.Node{idxNQ,1}{11,1}';
            children = Quadtree.getchildren(idxNQ);
            if ~isempty(children) % yes
                
                % select neighbour Quad children
                if dir == 1 % West NQ
                    idxC = children([3,4]); % NE and SE children
                    idxV = [2,3];
                    idxM = 4;
                elseif dir == 2 % South NQ
                    idxC = children([1,3]); % NW and NE children
                    idxV = [3,4];  
                    idxM = 1;
                elseif dir == 3 % East NQ
                    idxC = children([1,2]); % NW and SW children
                    idxV = [1,4];
                    idxM = 2;
                elseif dir == 4 % North NQ
                    idxC = children([2,4]); % SW and SE children
                    idxV = [2,1];
                    idxM = 3;
                end
                
                % get midPoint from neighbour Quad children
                pt1 = Quadtree.Node{idxC(1),1}{7,1}(:,idxV(1));
                pt2 = Quadtree.Node{idxC(2),1}{7,1}(:,idxV(2));
                
                % compare midPoints
                if norm(pt1-pt2) < 1e-10
                    numMidPoints = numMidPoints + 1;
                    midPoints(:,idxM) = pt1;
                    locMidPoints(idxM) = idxM;
                end 
                
            else % no
                continue
            end
            
        else % no
            continue
        end
        
    end
end

existMidPoints = midPoints(1,:) ~= -99;
midPoints = midPoints(:,existMidPoints);
locMidPoints = locMidPoints(existMidPoints);


end
