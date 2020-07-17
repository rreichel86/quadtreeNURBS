function s = findspan(n,p,u,U)
%Adapted from the NURBS Book


if abs(u-U(n+2))<1e-10, s=n; return,  end
% if abs(u-U(1))<1e-10;idx = find((U-u)>1e-10,1);s=idx-1;return;end 
% if abs(u-U(end))<1e-10;idx = find((U-u)>-1e-10,1);s=idx-1;return;end 
% idx = find((U-u)>1e-10,1);s=idx-1;
low = p;                                       
high = n + 1;                                   
mid = floor((low + high) / 2);                  
while (u < U(mid+1) || u >= U(mid+2))           
    if (u < U(mid+1))                           
        high = mid;                             
    else                                       
        low = mid;                               
    end  
    mid = floor((low + high) / 2);             
end                                            
                                               
s = mid;                                        
                                                
