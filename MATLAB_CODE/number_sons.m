function [nSons]= number_sons(Location,Loc_Current)
% number_sons function takes as input the address book of all the nodes and
% the position of the actual quad and checks how many levels of
% decompoition it has, by comparing its coordinates with the addresses of
% the book, which have the same beggining. Gives as output a scalar nSons
% beeing it the number of levels of decomposition

%finds elements with at least same length
idx = cellfun('length',Location)>=length(Loc_Current);
out = cell(length(Location),1);
for j=1:length(idx)
    if idx(j)==1
        %finds elements with the same depth as Loc_Current
        out{j}=Location{j}(1:length(Loc_Current));
    end
end
k=1;
for j=1:length(out)
    tf=isequal(out{j},Loc_Current);
    %finds the indexes of the descendents
    if tf;idx2(k)=j;
        k=k+1;
    end
end
Sons=cell(length(idx2),1);
%creates a cell with the descendents
for j=1:length(idx2);Sons{j}=Location{idx2(j)};end
%check how many desdendents
[s,d] = cellfun(@size,Sons);
nSons= max([d])-length(Loc_Current);
end