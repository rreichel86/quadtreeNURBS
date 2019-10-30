function [position_Q] = position(l,k,i,pos_aux,Quadtree)
% position function takes as input the geometrical definition of a quad and
% gives as ouput its position at [11 12 21] format

Location=Location_Quads(Quadtree);

if isempty(pos_aux)
if k == 1
    pos=[l];
end

if k>1
    if k>length(Location{end})
        pos=[Location{end} 1];
    end
    if k==length(Location{end})
        pos=[Location{end}];
        pos(end)=i;
    end
    if k<length(Location{end})
        pos=[Location{end}(1:k)];
        pos(end)=i;
    end
end
else
    pos=[pos_aux i];
end


for m=1:length(pos)
    
    if pos(m) == 1
        position_Q((2*(m-1)+1))=1; 
        position_Q((2*(m-1)+2))=1;
    end
    if pos(m) == 2
        position_Q((2*(m-1)+1))=2; 
        position_Q((2*(m-1)+2))=1;
    end
    if pos(m) == 3
        position_Q((2*(m-1)+1))=1; 
        position_Q((2*(m-1)+2))=2;
    end
    if pos(m) == 4
        position_Q((2*(m-1)+1))=2; 
        position_Q((2*(m-1)+2))=2;
    end

end
end