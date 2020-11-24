function locQ = ref2loc(refQ)
% ref2loc: compute Quad location in 1:4 format (11 = 1, 21 = 2, 12 = 3, 22 = 4)

lenlocQ = length(refQ)/2;
locQ = zeros(lenlocQ,1);

for l = 1:lenlocQ
    
    pos = refQ(2*l-1:2*l);
    
    if isequal(pos,[1 1]) % NW Quad
        loc = 1;
    elseif isequal(pos,[2 1]) % SW Quad
        loc = 2;
    elseif isequal(pos,[1 2]) % NE Quad
        loc = 3;
    elseif isequal(pos,[2 2]) % SW Quad
        loc = 4;
    end
    
    locQ(l) = loc;
    
end
end