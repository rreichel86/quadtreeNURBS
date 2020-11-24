function idxNQ = findNeighbour(Quadtree,idxQ,refQ,refNQ)

    lv = length(refNQ);
    lvNCA = 0;
    for l = lv:-2:2 
        if isequal( refQ(l-1:l), refNQ(l-1:l) )  
            lvNCA = l;
            break;  
        end        
    end    

    % nearest common ancesor idx and ref
    idxA = idxQ;
    for l = 1:(lv - lvNCA)/2
        idxA = Quadtree.Parent(idxA);     
    end

    % search neighbour idx 
    posNQ = refNQ(lvNCA+1:end);
    
    idxNQ = idxA;
    for i = 1:length(posNQ)/2
        
        pos = posNQ(2*i-1:2*i);
        
        if isequal(pos,[1 1]) % NW Quad
            loc = 1;
        elseif isequal(pos,[2 1]) % SW Quad
            loc = 2;
        elseif isequal(pos,[1 2]) % NE Quad
            loc = 3;
        elseif isequal(pos,[2 2]) % SW Quad
            loc = 4;
        end
        
        if idxNQ ~= 1
           children = Quadtree.Node{idxNQ,1}{11,1};
        else   
           children = Quadtree.Node{idxNQ,1}{2,1};
        end 
         
        if ~isempty(children)
            idxNQ = children(loc);
        else
            break;
        end
        
    end    
    
end 