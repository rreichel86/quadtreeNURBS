function ker = computePolygonKernel(poly,dbg)
% computePolygonKernel: compute kernel of given polygon
%
% INPUT:
% poly ----------------------- polygon's vertices given in CCW order
% poly = [x-coor, y-coor]
%             
% ker ------------------------ contains the vertices of the polygon's kernel
%                              only if the given polygon is star-shaped
%   
% -------------------------------------------------------------------------

if ~exist('dbg','var')
    dbg = 0;
end  

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

% plot polygon 
if dbg == 1
x = [ poly(:,1)', poly(1,1)'];
y = [ poly(:,2)', poly(1,2)'];  
plot(x,y, '-r')    
hold on 
end

n = size(poly,1); % number of vertices
% Identify concave vertices
zhl = 0;
cc = [];
for j = n:-1:1

    if j == 1 
        a = n;
    else
        a = j-1;
    end     
    
    if j == n
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
        cc(zhl,:) = j;
    end 

end 
ker = poly;

if ~isempty(cc)
    ncc = size(cc,1); % number of convace vertices
else     
    ker = poly;
    % plot polygon's kernel
    if dbg == 1
        x = [ poly(:,1)', poly(1,1)'];
        y = [ poly(:,2)', poly(1,2)'];
        patch(x,y,'green','FaceAlpha',0.5);
    end 
    return
end

% loop over convace vertices
for idx = 1:ncc 
    
    % A - vertex adjacent to current concave vertex (prev)
    % B - current concave vertex 
    % D - vertex adjacent to current concave vertex (next)
    if cc(idx) == 1
        Ai = poly(size(poly,1),:);
    else
        Ai = poly(cc(idx)-1,:);
    end
    
    Bi = poly(cc(idx),:);
    
    if cc(idx) == size(poly,1)
        Di = poly(1,:);
    else
        Di = poly(cc(idx)+1,:);
    end

    k = 1;
    while k <= ( size(ker,1) )
        
        % triangle A-B-C1
        % triangle A-B-C2
        A = Ai;
        B = Bi;
        C1 = ker(k,:);
        if k == size(ker,1)
            b = 1;
        else
            b = k+1;
        end
        C2 = ker(b,:);
        
        % compute orientation of A-B-C1 and A-B-C2 
        s1 = orientation(A,B,C1);
        s2 = orientation(A,B,C2);
        
        if s1 == 0 || s2 == 0
            k = k + 1;
            continue
        end
        
        % if the triangles A-B-C1 and A-B_C2 have opposite orientation
        % compute intersection U of the ray A-B and polygon edge C1-C2
        if ( s1 == -s2 )
            
            DD = [(B-A)' -(C2-C1)'];
            Dt = [(C1-A)' -(C2-C1)'];
            Ds = [(B-A)' (C1-A)'];
            
            tt = det(Dt)/det(DD);
            ss = det(Ds)/det(DD);
            
            tt*(B-A) + A;
            U = ss*(C2-C1) + C1;
            
            % insert U in Kernel between C1-C2
            if b == 1
                ker(size(ker,1)+1,:) = U;
            else
                temp = ker(b:end,:);
                ker(k+1,:) = U;
                ker(b+1:end+1,:) = temp;
            end
            
        end
        k = k + 1;
    end
    
    k = 1;
    while k <= ( size(ker,1) )
        
        % triangle B-D-C1
        % triangle B-D-C2
        A = Bi ;
        B = Di ;
        C1 = ker(k,:);
        if k == size(ker,1)
            b = 1;
        else
            b = k+1;
        end
        C2 = ker(b,:);
        
        % compute orientation of B-D-C1 and B-D-C2 
        s1 = orientation(A,B,C1);
        s2 = orientation(A,B,C2);
        
        if s1 == 0 || s2 == 0
            k = k + 1;
            continue
        end
        
        % if the triangles B-D-C1 and B-D_C2 have opposite orientation
        % compute intersection U of the ray B-D and polygon edge C1-C2
        if ( s1 == -s2)
            
            DD = [(B-A)' -(C2-C1)'];
            Dt = [(C1-A)' -(C2-C1)'];
            Ds = [(B-A)' (C1-A)'];
            
            tt = det(Dt)/det(DD);
            ss = det(Ds)/det(DD);
            
            tt*(B-A) + A;
            U = ss*(C2-C1) + C1;
            
            % insert U in Kernel between C1-C2
            if b == 1
                ker(size(ker,1)+1,:) = U;
            else
                temp = ker(b:end,:);
                ker(k+1,:) = U;
                ker(b+1:end+1,:) = temp;
            end
            
        end
        k = k + 1;
    end
    
    zhl = 0;
    tmp = [];
    for k = 1 : size(ker,1)
        C1 = ker(k,:);
        
        % compute orientation of A-B-C1 and B-D-C1
        s1 = orientation(Ai,Bi,C1);
        s2 = orientation(Bi,Di,C1);
        
        % if triangle A-B-C1 or triangle B-D-C1 is oriented CW
        % Store C1s
        if ( s1 == -1 || s2 == -1)
            
%             plot(C1(1),C1(2),'ro')
%             hold on
            
            zhl = zhl + 1;
            tmp(zhl) = k;
        end
    end
    
    % Delete C1s from the polygon's kernel
    ker(tmp,:)=[];
    if ~isempty(ker)
%         x = [ ker(:,1)', ker(1,1)'];
%         y = [ ker(:,2)', ker(1,2)'];
%         plot(x,y, '-')
%         hold on
    else
        break
    end
    
    
end
% end loop over concave vertices 

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
    % plot polygon's kernel
    if dbg == 1
        patch(x,y,'green','FaceAlpha',0.5);
    end 
end

end
