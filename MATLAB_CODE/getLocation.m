function [coor_nodes]=getLocation(xcoor1,xcoor2,ycoor1,ycoor2,grad)

ninode = grad - 1;
idx = [1:grad-1]';
coor_nodes = zeros(ninode,2);

derta_x = abs(xcoor1-xcoor2)/grad;
derta_y = abs(ycoor1-ycoor2)/grad;

    
if xcoor1 <= xcoor2
    coor_nodes(:,1) = xcoor1 + derta_x*idx;
else 
    coor_nodes(:,1) = xcoor1 - derta_x*idx;
end

if ycoor1 <= ycoor2
    coor_nodes(:,2) = ycoor1 + derta_y*idx;
else 
    coor_nodes(:,2) = ycoor1 - derta_y*idx;
end
    
                       
            
end
   