function [Quadtree] = QuadtreeMerge(Quadtree,seedingPoints)

% number of seeding points
nSeedingPoints = size(seedingPoints,1);
idxQs = [];

for i = 1:nSeedingPoints
    
    idx = 1;
    while true
        
%         if idx == 1
%             idxChildren = Quadtree.Node{idx,1}{2,1};
%         else
%             idxChildren = Quadtree.Node{idx,1}{11,1};
%         end

        idxChildren = Quadtree.getchildren(idx);
        
        if isempty(idxChildren)
            break
        end
        
        for j = 1:4
            
            Quad = Quadtree.Node{idxChildren(j),1}{10,1}(1:2,1:4);
            
            minCoords = [Quad(1,1),Quad(2,1)];
            maxCoords = [Quad(1,2),Quad(2,3)];
            
            ptInQuad = isPointInQuad(minCoords,maxCoords,seedingPoints(i,:));
            
            if ptInQuad == 1
                idx = idxChildren(j);
                break
            end
        end
        
    end
    idxQs = [idxQs idx];

end

idxQs = unique(idxQs);
nidxQs = length(idxQs);

i = nidxQs;

while  nidxQs ~= 0
    
    idx_Q = idxQs(i);
    Qsiblings = Quadtree.getsiblings(idx_Q);
    idx_Father = Quadtree.getparent(idx_Q);
    
    sib_members = ismember(Qsiblings,idxQs);
    count_sib_members = sum(sib_members);
    i = i - count_sib_members;
    
    if count_sib_members == 4
        
        data = Quadtree.get(idx_Father);
        data{11} = [];
        Quadtree = Quadtree.set(idx_Father,data);
        
        for j = 4:-1:1
            Quadtree = Quadtree.removenode(Qsiblings(j));
        end
        
        
    end     
    
    idxQs = setdiff(idxQs,Qsiblings);
    nidxQs = length(idxQs);
    
end    

end