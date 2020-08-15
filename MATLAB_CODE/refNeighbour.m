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