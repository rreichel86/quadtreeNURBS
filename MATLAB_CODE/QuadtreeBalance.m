function [Quadtree] = QuadtreeBalance(Quadtree, controlPoints,...
    knots, weights, degree,Boundary)

figure(1);


references = cellfun(@(Q) Q(2),Quadtree.Node);
references{1} = [];

% Quadtree leaves
idxLeaves = Quadtree.findleaves';
Leaves = Quadtree.Node(idxLeaves);
numLeaves = length(Leaves);

%for l = 1 : numLeaves
l = 1;
while l <= numLeaves
    
    idxLeaf = idxLeaves(l);
    idxLeaves(l) = 0;
    % current Quad
    % current Quad reference
    refLeaf = Leaves{l,1}{2,1}(1:end);
    % current Quad location in 1:4 format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
    locLeaf = ref2loc(refLeaf);
    leafSubdivs = countSubdivs(references,refLeaf);
    leafFather = Quadtree.Node{Quadtree.Parent(idxLeaf),1};
    idxLeafFather = Quadtree.Parent(idxLeaf);
    
    % plot current Quad
    LeafXcoor = Leaves{l,1}{10,1}(1,1:end);
    LeafYcoor = Leaves{l,1}{10,1}(2,1:end);
    patch(LeafXcoor', LeafYcoor', 'red','FaceAlpha',.2)
    hold on;
    
    % search for current Quad neighbour
    idxNeighbours = zeros(4,1);
    for dir = 1:4
        [exist_NQ, refNQ] = refNeighbour(refLeaf,dir);
        % neighbour Quad level
        levelNQ = length(refNQ);
        
        if exist_NQ == 1
            for level = levelNQ:-2:2
                idx = cellfun(@(x) isequal(x, refNQ(1:level)),references);
                
                if any(idx)
                    idxNQ = find(idx == 1);
                    idxNeighbours(dir) = idxNQ;
                    NQ = Quadtree.Node(idx);
                    break
                end
                
            end
            NQXcoor = NQ{1,1}{10,1}(1,1:end);
            NQYcoor = NQ{1,1}{10,1}(2,1:end);
            plot(NQXcoor', NQYcoor', 'b-')
            hold on;
            
            
            if splittQ(Quadtree,dir,idxNQ) == 1
                [Quadtree] = Decompose_balance(Quadtree,controlPoints, ...
                    knots, weights, degree,idxLeaf,locLeaf,idxLeafFather,Boundary);
                
                idxNewLeaves = Quadtree.Node{idxLeaf,1}{11,1}';
                newLeaves = Quadtree.Node(idxNewLeaves);
                idxLeaves = [idxLeaves; idxNewLeaves];
                Leaves(end+1:end+4) = newLeaves;
                
                numLeaves = numLeaves + 4;
                
                references = cellfun(@(Q) Q(2),Quadtree.Node);
                references{1} = [];
                
            end
        end
    end
    
    %
    
    
    
    
    l = l + 1;
    
end

end

function [exist_NQ, refNQ] = refNeighbour(refQ,dir)
% refNeighbour: determine reference of searched neighbour

posQ = refQ(end-1:end);
level = length(refQ)/2;
% convert into integer form N_x, N_y
N = zeros(level,2);
N(:,1) = refQ(2:2:end) - 1; % N_x
N(:,2) = refQ(1:2:end) - 1; % N_y
% compute lim_x, lim_y
lim = zeros(2,1);
lim(1) = calcLim(level, N(:,1)); % lim_x
lim(2) = calcLim(level, N(:,2)); % lim_y
% perform binary transformation to obtain the 4 possible edge neighbours
%   1 - West
%   2 - South
%   3 - East
%   4 - North
% and determine reference of searched edge neighbour by interweaving
% new N_x and N_y
[exist_NQ, refNQ] = edgeNeighbour(posQ, dir, level, N, lim);
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

splitt = 0;

children = Quadtree.Node{idxNQ,1}{11,1}';
if isempty(children)
    return
end

if dir == 1
    idx = children([3,4]);
elseif dir == 3
    idx = children([1,2]);
elseif dir == 2
    idx = children([1,3]);
elseif dir == 4
    idx = children([2,4]);
end

child_1 = Quadtree.Node{idx(1),1}{11,1};
child_2 = Quadtree.Node{idx(2),1}{11,1};
if ~isempty(child_1) || ~isempty(child_2)
    splitt = 1;
end

end
