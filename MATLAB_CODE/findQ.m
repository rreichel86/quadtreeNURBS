function idxQ = findQ(Quadtree,refQ)

   
    idxQ = 1;
    for i = 1:length(refQ)/2
        
        pos = refQ(2*i-1:2*i);
        
        if isequal(pos,[1 1]) % NW Quad
            loc = 1;
        elseif isequal(pos,[2 1]) % SW Quad
            loc = 2;
        elseif isequal(pos,[1 2]) % NE Quad
            loc = 3;
        elseif isequal(pos,[2 2]) % SW Quad
            loc = 4;
        end
        
        if idxQ ~= 1
           children = Quadtree.Node{idxQ,1}{11,1};
        else   
           children = Quadtree.Node{idxQ,1}{2,1};
        end 
         
        if ~isempty(children)
            idxQ = children(loc);
        else
            break;
        end
        
    end    
    
end 