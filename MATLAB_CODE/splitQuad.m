function [ children, vertices ] = splitQuad(quadVertices,aCCW)
% splitQuad: generate 4 children from the current parent Quad 

if ~exist('aCCW', 'var')
    aCCW = 0;
end     

if aCCW == 0
    
    xmin = min(quadVertices(1,:));
    xmax = max(quadVertices(1,:));
    
    ymin = min(quadVertices(2,:));
    ymax = max(quadVertices(2,:));
    
    % quad's vertices arranged CCW
    vertices = zeros(2,9);
    
    vertices(:,1) = [xmin; ymin];
    vertices(:,2) = [xmax; ymin];
    vertices(:,3) = [xmax; ymax];
    vertices(:,4) = [xmin; ymax];
else
    vertices = zeros(2,9);
    vertices(:,1:4) = quadVertices;
    
end 

% bi-linear shape functions
Nh = @(xsi,eta) shape_func_2d_lin(xsi,eta);

vertices(:,5) = vertices(:,1:4) * Nh(0,-1);
vertices(:,6) = vertices(:,1:4) * Nh(1,0);
vertices(:,7) = vertices(:,1:4) * Nh(0,1);
vertices(:,8) = vertices(:,1:4) * Nh(-1,0);
vertices(:,9) = vertices(:,1:4) * Nh(0,0);


children = [8,9,7,4;...
            1,5,9,8;...
            9,6,3,7;...
            5,2,6,9];

end

function Nh = shape_func_2d_lin(xsi,eta)

Nh = zeros(4,1);

Nh(1) = (1-xsi)*(1-eta)/4;
Nh(2) = (1+xsi)*(1-eta)/4;
Nh(3) = (1+xsi)*(1+eta)/4;
Nh(4) = (1-xsi)*(1+eta)/4;

end 
