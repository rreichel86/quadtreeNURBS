function [Quadtree] = QuadtreeBalance(Quadtree, controlPoints,...
    knots, weights, degree,Boundary)

figure(1);

references = cellfun(@(Q) Q(2),Quadtree.Node);
references{1} = [];

% array with Quadtree leaves indices
idxLeaves = Quadtree.findleaves';
numLeaves = length(idxLeaves);

l = 1;
while l <= numLeaves
    
    % get current Quad index
    idxLeaf = idxLeaves(l);
    idxLeaves(l) = 0;
    % current Quad
    % current Quad reference
    refLeaf = Quadtree.Node{idxLeaf,1}{2,1}(1:end);
    levelLeaf = length(refLeaf);
    % current Quad location in 1:4 format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
    locLeaf = ref2loc(refLeaf);
    % current Quad father index
    idxLeafFather = Quadtree.Parent(idxLeaf);
    
    % plot current Quad
    % LeafXcoor = Quadtree.Node{idxLeaf,1}{10,1}(1,1:end);
    % LeafYcoor = Quadtree.Node{idxLeaf,1}{10,1}(2,1:end);
    % patch(LeafXcoor', LeafYcoor', 'red','FaceAlpha',.2)
    % hold on;
    
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
            for level = levelNQ:-2:2
                idx = cellfun(@(x) isequal(x, refNQ(1:level)),references);
                if any(idx)
                    idxNQ = find(idx,1,'first');
                    % update array with found neighbours indices
                    niNQs = niNQs + 1;
                    iNQs(niNQs) = idxNQ;
                    break
                end
            end
            % plot found neighbour Quad
            % NQ = Quadtree.Node(idx);
            % NQXcoor = NQ{1,1}{10,1}(1,1:end);
            % NQYcoor = NQ{1,1}{10,1}(2,1:end);
            % plot(NQXcoor', NQYcoor', 'b-')
            % hold on;
            
            if is_splitted == 1
                continue
            end
            
            % check if current Quad has to be split
            if splittQ(Quadtree,dir,idxNQ) == 1
                % split current Quad 
                [Quadtree] = Decompose_balance(Quadtree,controlPoints, ...
                    knots, weights, degree,idxLeaf,locLeaf,idxLeafFather,Boundary);
                
                is_splitted = 1;
                idxNewLeaves = Quadtree.Node{idxLeaf,1}{11,1}';
                % insert current Quad children into array with leaves indices
                idxLeaves = [idxLeaves; idxNewLeaves];
                numLeaves = numLeaves + 4;
                
                references = cellfun(@(Q) Q(2),Quadtree.Node);
                references{1} = [];
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
            children = Quadtree.Node{iNQs(i),1}{11,1}';
            if isempty(children) % yes

                % Check if neighbour Quad is larger than 
                % the current Quad children 
                if (levelLeaf + 2) - levelNQ > 2
                    % insert neighbour Quad into array with leaves indices
                    idxLeaves = [idxLeaves; iNQs(i)];
                    numLeaves = numLeaves + 1;
                end
            end
        end
    end
    l = l + 1;
end
end



function lim = calcLim(level, N)
% calcLim: calculate lim
% working backwards through N, this is the level
% at which N(i) first becomes not equal to N(i-1).

for l = level:-1:1
    if (l == 1)
        lim = l;
        break;
    elseif N(l) ~= N(l-1)
        lim = l;
        break;
    end
end
end

function [exist_NQ, refNQ] = edgeNeighbour(posQ, dir, level, N, lim)
% edgeNeighbour: perform binary transformation to obtain the
% 4 possible edge neighbours:
%   1 - West
%   2 - South
%   3 - East
%   4 - North
% and determine reference of searched edge neighbour by interweaving
% new N_x and N_y

has_same_father = 1;

% Check if the neighbour (dir) we are looking for has the same father
% as the current Quad
if isequal(posQ,[1 1]) % NW Quad
    if dir == 1 || dir == 4
        has_same_father = 0;
    end
elseif isequal(posQ,[2 1]) % SW Quad
    if dir == 1 || dir == 2
        has_same_father = 0;
    end
elseif isequal(posQ,[1 2]) % NE Quad
    if dir == 3 || dir == 4
        has_same_father = 0;
    end
elseif isequal(posQ,[2 2]) % SW Quad
    if dir == 3 || dir == 2
        has_same_father = 0;
    end
end
% has the same father ?
if has_same_father == 1 % yes
    % peform binary operation on N(level)
    if dir == 1 || dir == 3
        [exist_NQ, N(:,1)] = binaryTransformation(level, N(:,1));
    elseif dir == 2 || dir == 4
        [exist_NQ, N(:,2)] = binaryTransformation(level, N(:,2));
    end
else % no
    % perform binary operation on N(level) to N(lim-1)
    if dir == 1 || dir == 3
        [exist_NQ, N(:,1)] = binaryTransformation(level, N(:,1), lim(1));
    elseif dir == 2 || dir == 4
        [exist_NQ, N(:,2)] = binaryTransformation(level, N(:,2), lim(2));
    end
end

% determine reference of searched neighbour by interweaving new N_x and N_y
if exist_NQ == 1
    refNQ = zeros(1,2*level);
    refNQ(2:2:end) = N(:,1) + 1; % N_x
    refNQ(1:2:end) = N(:,2) + 1; % N_y
else
    refNQ = -1;
end
end

function [status, N] = binaryTransformation(level, N, lim)
% binaryTransformation: perform binary operation on N

status = 1;

if ~exist('lim','var')
    N(level) = 1 - N(level); % on N(level)
else
    if lim == 1
        status = 0;
        return
    end
    N(level:-1:lim-1) = - N(level:-1:lim-1) + 1; % on N(level) to N(lim-1)
end
end


function locQ = ref2loc(refQ)

lenlocQ = length(refQ)/2;
for l = 1:lenlocQ
    
    pos = refQ(2*l-1:2*l);
    
    if isequal(pos,[1 1]) % NW Quad
        loc = 1;
    elseif isequal(pos,[2 1]) % SW Quad
        loc = 2;
    elseif isequal(pos,[1 2]) % NE Quad
        loc = 3;
    elseif isequal(pos,[2 2]) % SW Quad
        loc = 4;
    end
    
    locQ(l) = loc;
    
end
end


function splitt = splittQ(Quadtree,dir,idxNQ)
% splittQ: check if Quad has to be splitted

splitt = 0;

% has neighbour Quad children ?
children = Quadtree.Node{idxNQ,1}{11,1}';
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
child_1 = Quadtree.Node{idx(1),1}{11,1};
child_2 = Quadtree.Node{idx(2),1}{11,1};
if ~isempty(child_1) || ~isempty(child_2)
    splitt = 1;
end

end
