function [Quadtree] = QuadtreeBalance(Quadtree,NURBS)
% QuadtreeBalance: balance quadtree, apply 2:1 rule
%
% INPUT:
% Quadtree ------------------- Quadtree data structure
% NURBS definition
% NURBS.degree --------------- NURBS degree
% NURBS.knots ---------------- NURBS knot vector
% NURBS.controlPoints -------- NURBS control points
% NURBS.weights -------------- NURBS weights
% 
% OUTPUT:
% Quadtree ------------------- Quadtree data structure
% 
% -------------------------------------------------------------------------

% array with Quadtree leaves indices
idxLeaves = Quadtree.findleaves';

l = 1;
while l <= length(idxLeaves)
    
    % get current Quad index
    idxLeaf = idxLeaves(l);
    idxLeaves(l) = 0;
    % current Quad
    % current Quad reference
    refLeaf = Quadtree.Node{idxLeaf,1}{2,1}(1:end);
    levelLeaf = length(refLeaf);
    % current Quad location in 1:4 format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
    % locLeaf = ref2loc(refLeaf);
    % current Quad father index
    % idxLeafFather = Quadtree.Parent(idxLeaf);
    
    is_splitted = 0;
    % array with found neighbours indices
    iNQs = zeros(4,1);
    % count found neighbours
    niNQs = 0;
    % search for current Quad neighbours
    % Loop over directions:
    % 1 - West
    % 2 - South
    % 3 - East
    % 4 - North
    for dir = 1:4
        % determine possible neighbour Quad in current direction 
        [exist_NQ, refNQ] = refNeighbour(refLeaf,dir);
        if exist_NQ == 1
            % look for neighbour Quad reference in reference array
            % and get its position or its ancestor position
            levelNQ = length(refNQ);
            
            idxNQ = findNeighbour(Quadtree,idxLeaf,refLeaf, refNQ);
            % update array with found neighbours indices
            niNQs = niNQs + 1;
            iNQs(niNQs) = idxNQ;
                        
            if is_splitted == 1
                continue
            end
            
            % check if current Quad has to be split
            if splittQ(Quadtree,idxLeaf,dir,idxNQ) == 1
                % split current Quad 
                [Quadtree] = Decompose_helper(Quadtree,NURBS,idxLeaf);
                
                is_splitted = 1;
                % idxNewLeaves = Quadtree.Node{idxLeaf,1}{11,1}';
                idxNewLeaves = Quadtree.getchildren(idxLeaf)';
                % insert current Quad children into array with leaves indices
                idxLeaves = [idxLeaves; idxNewLeaves];
                
            end
        end
    end

    % if current quad was splitted
    % check if current Quad has neighbours that need to be split.
    if is_splitted == 1
        
        % Loop over array with indices of found neighbours
        for i = 1: niNQs
            
            refNQ = Quadtree.Node{iNQs(i),1}{2,1}(1:end);
            % neighbour Quad level
            levelNQ = length(refNQ);
            % is neighbour Quad a leaf?
            % children = Quadtree.Node{iNQs(i),1}{11,1}';
            children = Quadtree.getchildren(iNQs(i));
            if isempty(children) % yes

                % Check if neighbour Quad is larger than 
                % the current Quad children 
                if (levelLeaf + 2) - levelNQ > 2
                    % insert neighbour Quad into array with leaves indices
                    idxLeaves = [idxLeaves; iNQs(i)];

                end
            end
        end
    end
    l = l + 1;
end
end


function splitt = splittQ(Quadtree,idxQ,dir,idxNQ)
% splittQ: check if Quad has to be splitted

splitt = 0;
% check if current Quad is a leaf
% children = Quadtree.Node{idxQ,1}{11,1}';
children = Quadtree.getchildren(idxQ);
if ~isempty(children)
    return % no
end

% has neighbour Quad children ?
% children = Quadtree.Node{idxNQ,1}{11,1}';
children = Quadtree.getchildren(idxNQ);
if isempty(children)
    return % no
end
% neighbour Quad child children to be checked
if dir == 1 % West NQ
    idx = children([3,4]); % NE and SE children
elseif dir == 3 % East NQ
    idx = children([1,2]); % NW and SW children
elseif dir == 2 % South NQ
    idx = children([1,3]); % NW and NE children
elseif dir == 4 % North NQ
    idx = children([2,4]); % SW and SE children
end
% check if the corresponding neighbour Quad child has also children
% child_1 = Quadtree.Node{idx(1),1}{11,1};
child_1 = Quadtree.getchildren(idx(1));
% child_2 = Quadtree.Node{idx(2),1}{11,1};
child_2 = Quadtree.getchildren(idx(2));
if ~isempty(child_1) || ~isempty(child_2)
    splitt = 1;
end

end
