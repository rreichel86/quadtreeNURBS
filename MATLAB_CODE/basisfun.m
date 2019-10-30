function N = basisfun(i,u,p,U)                 
%Adapted from the NURBS Book

i = i + 1; 
                                                 
left = zeros(p+1,1);                              
right = zeros(p+1,1);                            
                                                
N(1) = 1;                                         
for j=1:p                                         
    left(j+1) = u - U(i+1-j);                     
    right(j+1) = U(i+j) - u;                      
    saved = 0;                                    
    for r=0:j-1                                   
        temp = N(r+1)/(right(r+2) + left(j-r+1));  
        N(r+1) = saved + right(r+2)*temp;         
        saved = left(j-r+1)*temp;                 
    end                                            
    N(j+1) = saved;                               
end                                               
   
                                                  
