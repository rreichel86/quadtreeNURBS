function ker = computePolygonKernel(poly)
% clear all
%    p = [20.4 13.5;...
%         5.6 40.9;...
%         24.7 49.2;...
%         30.2 69.5;...
%         67.4 64.3;...
%         59.9 48.4;...
%         63.6 27.6;...
%         43.5 31.6;...
%         41 18.8];
%       
    
%     p = [14.4 7.2;...
%          9.2 22.2;...
%         23.2 32.1;...
%         14.7 48.4;...
%         30.8 64.4;...
%         59.2 55.9;...
%         53.9 33.3;...
%         73.7 22.3;...
%         53.1 12.5;...
%         41.8 21.8];
    
%      p = [41.2 6;...
%          18.6 9.3;...
%         20.4 31.9;...
%         7.1 21.8;...
%         4.1 38.4;...
%         17.4 56.3;...
%         38.3 51;...
%         61.9 62.3;...
%         61.9 22.9;...
%         44.8 31.9];
    
%     p = [7.7 39.9;...
%          17 59.7;...
%         30.8 43.2;...
%         51.2 62;...
%         65.3 41.6;...
%         48.7 43.7;...
%         62 26;...
%         41.4 10;...
%         22.6 26.3;...
%         15 6.5];
%     
% n = size(p,1);    
% for  i = 1 : n
%     
%     poly(n-i+1,:) = p(i,:);
%     
% end 


% area = poly(n,1) * poly(1,2) - poly(1,1) * poly(n,2);
% 
% for i = 1: n-1
% 
%  area = area + poly(i,1) * poly(i+1,2);
%  area = area - poly(i+1,1) * poly(i,2);
% 
% end 
% 
% if area > 0   
%     OP = 1;
% else 
%     OP = -1; 
% end

n = size(poly,1);

zhl = 0;
cc = [];
for j = n:-1:1


    if (j-1) == 0 
        a = n;
    else
        a = j-1;
    end     
    
    if (j+1) == n+1
        b = 1;
    else
        b = j+1;
    end 
    
    A0 = poly(a,:);
    B0 = poly(j,:);
    C0 = poly(b,:);

    s = orientation(A0,B0,C0);
    
    if s == -1 
        zhl = zhl + 1;
        cc(zhl,:) =j;
    end 

end 
ker = poly;

if ~isempty(cc)
    ncc = size(cc,1);
    
else     
    ker = poly;
    return
end

x = [ poly(:,1)', poly(1,1)'];
y = [ poly(:,2)', poly(1,2)'];  


% plot(x,y, '-r')    

% hold on 




for idx = 1:ncc 
    
cc_i = poly(cc(idx),:);    

    if cc(idx) == 1
        Ai = poly(size(poly,1),:);
    else
        Ai = poly(cc(idx)-1,:);
    end
    
    Bi = poly(cc(idx),:);
    
    if cc(idx) >= size(poly,1)
        Di = poly(1,:);
    else
        Di = poly(cc(idx)+1,:);
    end
    

for k = 1:size(ker,1)            
  
    %v(i-1)v(i)s(j)
    %v(i-1)v(i)s(j+1)
    A = Ai ;
    B = Bi ;
    C1 = ker(k,:);
    if (k+1) > size(ker,1)
        b = 1;
    else
        b = k+1;
    end 
       C2 = ker(b,:);
         
    s1 = orientation(A,B,C1);
    s2 = orientation(A,B,C2);

    if ( s1 ~= 0 && s2 ~= 0 && s1 == -s2) 
        k;
        b;
        DD = [(B-A)' -(C2-C1)'];
        Dt = [(C1-A)' -(C2-C1)'];
        Ds = [(B-A)' (C1-A)'];
        
        tt = det(Dt)/det(DD);
        ss = det(Ds)/det(DD);
        
        tt*(B-A) + A;
        U = ss*(C2-C1) + C1;
        
       
        if b == 1 
            ker(size(ker,1)+1,:) = U;
        else 
            temp = ker(b:end,:);
            ker(k+1,:) = U;
            ker(b+1:end+1,:) = temp;
        end 

    end     
end


for k = 1:size(ker,1)      
  
    %v(i)v(i+1)s(j)
    %v(i)v(i+1)s(j+1)
    A = Bi ;
    B = Di ;
    C1 =ker(k,:);
    if (k+1) > size(ker,1)
        b = 1;
    else
        b = k+1;
    end 
       C2 = ker(b,:);
       
    s1 = orientation(A,B,C1);
    s2 = orientation(A,B,C2);

    if ( s1 ~= 0 && s2 ~= 0 && s1 == -s2) 
        k;
        b;
        DD = [(B-A)' -(C2-C1)'];
        Dt = [(C1-A)' -(C2-C1)'];
        Ds = [(B-A)' (C1-A)'];
        
        tt = det(Dt)/det(DD);
        ss = det(Ds)/det(DD);
        
        tt*(B-A) + A;
        U = ss*(C2-C1) + C1;
        
       
        if b == 1 
            ker(size(ker,1)+1,:) = U;
        else 
            temp = ker(b:end,:);
            ker(k+1,:) = U;
            ker(b+1:end+1,:) = temp;
        end 

    end     
    
    
end


zhl = 0;
tmp = [];
for k = 1 : size(ker,1)
   C1 = ker(k,:);
                 
   s1 = orientation(Ai,Bi,C1);
   s2 = orientation(Bi,Di,C1);
          
   if ( s1 == -1 || s2 == -1)
       zhl = zhl + 1;
       k;      
       tmp(zhl) = k;
   end    
end

ker(tmp,:)=[];
if ~isempty(ker)
    x = [ ker(:,1)', ker(1,1)'];
    y = [ ker(:,2)', ker(1,2)'];  

%     plot(x,y, '-')    
%     hold on
else
    break
end     


end
ker;
nvertices = size(ker,1);

% check if kernel is not empty or 
% consist of at least 3 noncollinear vertices
if nvertices < 3
    ker = [];
elseif nvertices == 3
    if orientation(ker(1,:),ker(2,:),ker(3,:)) == 0
        ker = [];
    end
end

if ~isempty(ker)
   %patch(x,y,'green');
end

end
