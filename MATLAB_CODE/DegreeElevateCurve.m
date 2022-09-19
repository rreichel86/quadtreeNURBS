function [knots_new, Q, weights] = DegreeElevateCurve(degree,knots,controlPoints,weights,t)

% raise the degree from d to d+t, t>=1.
% Adapted from the NURBS book (Algorithm A5.9)



[mc,nc] = size(controlPoints);

n = nc - 1;

d = degree;
bezalfs =  zeros(d+1,d+t+1);   % coefficients for degree elevating the Bezier segments;                        
bpts = zeros(mc,d+1);          % dth-degree Bezier control points of the current segment;                           
ebpts = zeros(mc,d+t+1);       % (d + t)th-degree Bezier control points of the current segment;                           
Nextbpts = zeros(mc,d+1);      % leftmost control points of the next Bezier segment;                         
alfs = zeros(d,1);             % knot insertion alpha-times.

m = n+d+1;
ph = d+t;                                                                                                                                               
ph2 = floor(ph / 2); 

% Control points before order elevation
Pw = [controlPoints(1,:) .* weights ;controlPoints(2,:) .* weights; weights];
% knot vector after order eleveation
reKnots0 = sort(knots);
reKnots1 = unique(reKnots0);
numKnots_new = length(reKnots0) + t*length(reKnots1);
knots_new = zeros(1,numKnots_new);
% Control points after order elevation
nq = numKnots_new-(d+t)-1;
Qw = zeros(mc+1,nq);
                                                          
%compute bezier degree elevation coefficeients                                                          
bezalfs(1,1) = 1;                                         
bezalfs(d+1,ph+1) = 1;                                   

for i=1:ph2                                               
   inv = 1/bincoeff(ph,i);                                
   mpi = min(d,i);                                                                                               
   for j=max(0,i-t):mpi                                   
       bezalfs(j+1,i+1) = inv*bincoeff(d,j)*bincoeff(t,i-j); 
   end                                                       
end                                                      
                                                          
for i=ph2+1:ph-1                                          
   mpi = min(d,i);                                        
   for j=max(0,i-t):mpi                                   
       bezalfs(j+1,i+1) = bezalfs(d-j+1,ph-i+1);          
   end                                                       
end                                                       
                                                          
mh = ph;                                                                                                    
ua = knots(1);                                                
                                                          
for ii=0:mc                                            
   Qw(ii+1,1) = Pw(ii+1,1);                                
end                                                       
for i=0:ph                                                
   knots_new(i+1) = ua;                                          
end                                                       

kind = ph+1; 
cind = 1; 
                                                  

%initialise first bezier seg
lbz = 1;
r = -1; 
for i=0:d                                                 
   for ii=0:mc                                      
      bpts(ii+1,i+1) = Pw(ii+1,i+1);                      
   end                                                       
end 

a = d;                                                    
b = d+1;  

%big loop thru knot vector                                                         
while b < m                                               
   i = b;                                                 
   while b < m && knots(b+1) == knots(b+2)                        
      b = b + 1;                                          
   end                                                    
   mul = b - i + 1;                                       
   mh = mh + mul + t;                                    
   ub = knots(b+1);                                           
   oldr = r;                                              
   r = d - mul;                                           
                                                          
   %insert knot u(b) r times                                                       
   if oldr > 0                                            
      lbz = floor((oldr+2)/2);                           
   else                                                  
      lbz = 1;                                            
   end                                                    
   
   if r > 0                                              
      rbz = ph - floor((r+1)/2);                          
   else                                                   
      rbz = ph;                                           
   end                                                    
   
   if r > 0                                               
      %insert knot to get bezier segment                                                    
      numer = ub - ua;                                    
      for q=d:-1:mul+1                                    
         alfs(q-mul) = numer / (knots(a+q+1)-ua);            
      end                                           
      
      for j=1:r                                           
         save = r - j;                                    
         s = mul + j;                                     
                                                          
         for q=d:-1:s                                     
            for ii=0:mc                              
               tmp1 = alfs(q-s+1)*bpts(ii+1,q+1); 
               tmp2 = (1-alfs(q-s+1))*bpts(ii+1,q); 
               bpts(ii+1,q+1) = tmp1 + tmp2;           
            end                                              
         end                                              
         
         for ii=0:mc                                
            Nextbpts(ii+1,save+1) = bpts(ii+1,d+1);      
         end                                                 
      end                                                
   end   %end of insert knot                                                 
                                                          
                                                          
      
   if r > 0 
       RBZ = rbz+1;
   else
       RBZ = rbz;
   end

   %degree elevate bezier
   for i=lbz:RBZ   % Only points lbz,....,RBZ are used below                                      
      for ii=0:mc                                      
         ebpts(ii+1,i+1) = 0;                             
      end                                                    
      mpi = min(d, i);                                   
      for j=max(0,i-t):mpi                               
         for ii=0:mc                                  
            tmp1 = ebpts(ii+1,i+1); 
            tmp2 = bezalfs(j+1,i+1)*bpts(ii+1,j+1);
            ebpts(ii+1,i+1) = tmp1 + tmp2;                
         end                                                 
      end                                                    
   end   %end of degree elevating bezier                                                   
                                                          
                                                          
   if oldr > 1  %must remove knot u=k[a] oldr times  

      first = kind - 2;                                                  
      last = kind;                                       
      den = ub - ua;                                     
      bet = floor((ub-knots_new(kind)) / den);                   
                                                          
                                                          
      for tr=1:oldr-1    %knot removal loop                                 
         i = first;                                      
         j = last;                                        
         kj = j - kind + 1;                               
         while j-i > tr  %loop and compute the new control points 
                         %for one removal step
                                                          
                                                         
            if i < cind                                   
               alf = (ua-knots_new(i+1))/(ub-knots_new(i+1));          
               for ii=0:mc                             
                  tmp1 = alf*Qw(ii+1,i+1); 
                  tmp2 = (1-alf)*Qw(ii+1,i); 
                  Qw(ii+1,i+1) = tmp1 + tmp2;            
               end                                           
            end                                          
            if j >= lbz                                   
               if j-tr <= kind-ph+oldr                   
                  gam = (ub-knots_new(j-tr+1)) / den;            
                  for ii=0:mc                          
                     tmp1 = gam*ebpts(ii+1,kj+1); 
                     tmp2 = (1-gam)*ebpts(ii+1,kj+2); 
                     ebpts(ii+1,kj+1) = tmp1 + tmp2;     
                  end                                     
               else                                      
                  for ii=0:mc                         
                     tmp1 = bet*ebpts(ii+1,kj+1);                                     
                     tmp2 = (1-bet)*ebpts(ii+1,kj+2);                                     
                     ebpts(ii+1,kj+1) = tmp1 + tmp2;     
                  end                                        
               end                                        
            end                                          
            i = i + 1;                                    
            j = j - 1;                                    
            kj = kj - 1;                                  
         end                                            
                                                          
         first = first - 1;                               
         last = last + 1;                                 
      end                                                
   end   %end of removing knot n=k[a]                                                
                                                          
                                                      
   %load the knot ua                                                       
   if a ~= d                                            
      for i=0:ph-oldr-1                                   
         knots_new(kind+1) = ua;                                 
         kind = kind + 1;                                 
      end
   end                                                   
                                                         
   %load ctrl pts into ic                                                    
   for j=lbz:rbz 
       for ii=0:mc
           Qw(ii+1,cind+1) = ebpts(ii+1,j+1);            
       end                                                 
       cind = cind + 1;                                
   end  

   %setup for next pass thru loop                                                  
   if b < m
       for j=0:r-1 
           for ii=0:mc
               bpts(ii+1,j+1) = Nextbpts(ii+1,j+1);       
           end                                           
       end                                              
       for j=r:d                                        
           for ii=0:mc                              
               bpts(ii+1,j+1) = Pw(ii+1,b-d+j+1);         
           end                                              
       end                                                  
       a = b;                                           
       b = b+1;                                         
       ua = ub;                                         
                                                    
   else                                                       
       for i=0:ph 
           knots_new(kind+i+1) = ub;                           
       end                                                 
   end                                                    
end  % End big while loop       

Q = Qw(1:2,:) ./ Qw(3,:);
weights = Qw(3,:);

end 
                                                          
                                         

function b = bincoeff(n,k)
%  Computes the binomial coefficient.
%
%      ( n )      n!
%      (   ) = --------
%      ( k )   k!(n-k)!
%
%  b = bincoeff(n,k)
%
%  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.

                                                          
b = floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));     
end                                                      

function f = factln(n)
% computes ln(n!)
if n <= 1, f = 0; return, end
f = gammaln(n+1); %log(factorial(n));
end

