function f = lfunc(px,py,x1,y1,x2,y2)
    
    
    ax = (x2-x1);
    ay = (y2-y1);
    
    bx = (px - x1);
    by = (py - y1);
    
    f = -((ax*by - ay*bx)/(sqrt(ax*ax+ay*ay)));
end
    
 
  





