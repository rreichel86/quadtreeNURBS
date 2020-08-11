function [nSons]= number_sons(refList,currentRef)
% number_sons function takes as input the address book of all the nodes and
% the position of the actual quad and checks how many levels of
% decompoition it has, by comparing its coordinates with the addresses of
% the book, which have the same beggining. Gives as output a scalar nSons
% beeing it the number of levels of decomposition

% current Quad depth
depth = length(currentRef);
% get Quad, that are at least at the same depth
idx = cellfun('length',refList) >= depth;
filterRefList = refList(idx);
idx = cellfun(@(x) isequal(x(1:depth),currentRef),filterRefList);
% list descendents 
descendents = filterRefList(idx);
% check how many descendents
[~,d] = cellfun(@size,descendents);
nSons = max(d)-depth;
end
