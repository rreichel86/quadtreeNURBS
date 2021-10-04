function [Quadtree] = QuadtreeSplit(Quadtree,NURBS,seedingPoints)

% number of seeding points
nSeedingPoints = size(seedingPoints,1);
splitted = [];
fatherWasSplitted =  0;

for i = 1:nSeedingPoints
    
    idx = 1;
    fatherWasSplitted =  0;
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
  
    idxFather = Quadtree.Parent(idx);
    
    if  any (splitted == idxFather)
        fatherWasSplitted =  1;
    end
    
    if fatherWasSplitted ~= 1
        [Quadtree] = Decompose_helper(Quadtree,NURBS,idx);
        splitted = [splitted idx];
    end
    
end


end