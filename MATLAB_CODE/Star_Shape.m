function [Quadtree]=Star_Shape(Quadtree,controlPoints,NURBS, knots, weights, degree, Boundary)

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
l=Quadtree.findleaves();
flag = 0;


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
                [Quadtree]=starAlgorithm1(Quadtree,l(i));
                
                % Checks if kernels have already been found for that leaf
                % and breaks if true
                if ~isempty(Quadtree.Node{l(i),1}{12,1}) && ~isempty(Quadtree.Node{l(i),1}{13,1})
                    break
                end
                % Calls Algorithm2 for kerenl searching if the loop wasn't
                % broken
                [Quadtree]=starAlgorithm2(Quadtree,l(i));
                
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
    
end

end
