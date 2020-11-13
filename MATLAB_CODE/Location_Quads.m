function [Location]=Location_Quads(Quadtree)
% Location_Quads function takes as input the Quadtree and creates an
% address book in 1:4 format
l=length(Quadtree.Node);
% Location = cell(l,1);

% for i = 2:l
%     Current = Quadtree.Node{i,1}{2,1}(1:end);
%     for j = 1:length(Current)/2
%         aux = Current(2*j-1:2*j);
%         if isequal(aux,[1 1])
%             Location{i}=[Location{i} 1];
%         elseif isequal(aux,[2 1])
%             Location{i}=[Location{i} 2];
%         elseif isequal(aux,[1 2])
%             Location{i}=[Location{i} 3];
%         elseif isequal(aux,[2 2])
%             Location{i}=[Location{i} 4];
%         end                  
%     end    
% end

Location = cell(1,1);
Current = Quadtree.Node{l,1}{2,1}(1:end);
for j = 1:length(Current)/2
    aux = Current(2*j-1:2*j);
    if isequal(aux,[1 1])
        Location{1}=[Location{1} 1];
    elseif isequal(aux,[2 1])
        Location{1}=[Location{1} 2];
    elseif isequal(aux,[1 2])
        Location{1}=[Location{1} 3];
    elseif isequal(aux,[2 2])
        Location{1}=[Location{1} 4];
    end
end



end