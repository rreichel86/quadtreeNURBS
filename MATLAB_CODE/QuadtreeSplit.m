function [Quadtree] = QuadtreeSplit(Quadtree,NURBS,seedingPoints)

% number of seeding points
nSeedingPoints = size(seedingPoints,1);

for i = 1:nSeedingPoints
    
    idx = 1;
    while true
        
        idxChildren = Quadtree.getchildren(idx);
        
        if isempty(idxChildren)
            break
        end
        
        for j = 1:4
            
            Quad = Quadtree.Node{idxChildren(j),1}{7,1}(1:2,1:4);
                        
            if isPointInPolygon(Quad', seedingPoints(i,:))  ~= -1
                 
                idx = idxChildren(j);
                break
            end
        end
        
    end
    
    %     idxRef = Quadtree.Node{idx,1}{2,1};
    %     idxFather = Quadtree.Parent(idx);
    %     idxLoc = ref2loc(idxRef);
    
    [Quadtree] = Decompose_helper(Quadtree,NURBS,idx);
    
end


end