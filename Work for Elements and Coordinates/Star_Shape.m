function [Quadtree,element,nel]=Star_Shape(Quadtree,controlPoints,NURBS, knots, weights, degree, Boundary)

% Star_Shape is the main function which controls the kernel searching,
% First, it calls starAlgorithm1 for all leaves, that have a piece of the
% NURBS curve. If nessecery, starAlgorithm2 is called. If starAlgorithm2
% also fails to find kernels for both subdomains, a decompositon is performed 
% on the quad, for whcih no kernels can be found is performed

% Input: 
% Quadtree 
% data
% controlPoints 
% NURBS
% knots
% weights
% degree
% Boundary

% Output: 
% Quadtree: Updated Quadtree, containing the coordinates of vertices of the
% kernels at rows 12 and 13 for each leaf node in Quadtree

%%%%%% DELETE
%close all
%%%%%%%%%%%

% Find the leaves
l=Quadtree.findleaves()
flag = 0;

element_1 = cell(length(l),1);
% Loop through all leaves in searching for kernels
for i=1:length(l)
    
    % Obtains the interesections of the NURBS in the quad with the quad 
    intersections=Quadtree.Node{l(i),1}{5,1};
    % If empty or just 1 intersection, it will skip and won't search for
    % kernel in that quad
    if ~isempty(intersections) && length(intersections) ~= 1
        
        while true
            intersections=Quadtree.Node{l(i),1}{5,1};
            
            if ~isempty(intersections) && length(intersections) ~= 1
                % Calls Algorithm1 for kerenl searching
                [Quadtree,polygon1, polygon2]=starAlgorithm1(Quadtree,l(i));
                
                % Checks if kernels have already been found for that leaf
                % and breaks if true
                if ~isempty(Quadtree.Node{l(i),1}{12,1}) && ~isempty(Quadtree.Node{l(i),1}{13,1})
                    break
                end
                % Calls Algorithm2 for kerenl searching if the loop wasn't
                % broken
                [Quadtree,polygon1, polygon2]=starAlgorithm2(Quadtree,l(i));
                
                % Checks if kernels have already been found for that leaf
                % and breaks if true
                if ~isempty(Quadtree.Node{l(i),1}{12,1}) && ~isempty(Quadtree.Node{l(i),1}{13,1})
                    break
                end
                
                % If the loop still hasn't been broken, no kerenls were
                % found after Algorithm2 , so a decomposition is needed
                if true
                    [Quadtree]=Decompose_Star(Quadtree,NURBS,controlPoints, ...
                        knots, weights, degree,l(i),Boundary);
                    flag = 1;
                    break
                end
                
            else
                
                break
                
            end
            
        end
        
        if flag == 1
            % If the quad has been decemposed, the kerenl searching must
            % start from the beginning for the new situation 
            Star_Shape(Quadtree,controlPoints,NURBS, knots, weights, degree, Boundary);
            break;
            
        end
    end
    
    if isempty(intersections) || length(intersections) == 1
        
        quad=Quadtree.Node{l(i),1}{10,1}(1:2,1:4)
        
        element = [quad];
        
    elseif length(polygon1)>length(polygon2)
            element=[polygon1;polygon2,zeros(2,length(polygon1)-length(polygon2))]
    elseif length(polygon2)>length(polygon1)
            element=[polygon1,zeros(2,length(polygon2)-length(polygon1));polygon2]
        else
            element=[polygon1;polygon2]
        end
        
       element_1{i} = element; 
    end

nel=0

for i=1:length(l);
    intersections=Quadtree.Node{l(i),1}{5,1};
    
    if isempty(intersections) || length(intersections) == 1
        nel=nel+1
    else
        nel=nel+2
    end
end


count=1;
element=zeros(2*nel,1);

for i=1:length(element_1)
     [at,bt]=size(element_1{i});
 
    aut=count+at;
        element(count:aut-1,1:bt)=element_1{i};
        count=aut;
   
end

 end
