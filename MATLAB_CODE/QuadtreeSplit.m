function [Quadtree] = QuadtreeSplit(Quadtree,NURBS,seedingPoints)

% number of seeding points
nSeedingPoints = size(seedingPoints,1);

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
            
            Quad = Quadtree.Node{idxChildren(j),1}{7,1}(1:2,1:4);
            
            minCoords = [Quad(1,1),Quad(2,1)];
            maxCoords = [Quad(1,2),Quad(2,3)];
            
            ptInQuad = isPointInQuad(minCoords,maxCoords,seedingPoints(i,:));
            
            if ptInQuad == 1
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