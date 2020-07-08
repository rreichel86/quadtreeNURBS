function [polygon1,polygon2] = starSplitPolygons(Quad_current,points,Int)

% starSplitPolygons splits a domain, divided by a control polygon or
% polygonised NURBS curve into two subdomains

points_flip = fliplr(points); % Flip the values of points


% Stores the 4 vertices of the quad in variable polygon, in which later the
% intersection points will be inserted, preserving a clockwise direction of
% rotation, starting from the bottom left vertex
polygon = Quad_current; 
% Initialize variables to be used in the for loop
insertAfterPos = [];
insertions = 0;

% Code to insert intersection points in quad polygon
% loop trough the intersection points (exactly 2)
for l = 1:2
    flag1 = 0; % for controling when to brake the loop
    flag2 = 0; % for controling when to brake the loop
    % loop trough the lenght of polygon
    for i = 1:length(polygon)
        % Check if x coordinates of l-th intersection point and i-th vertex
        % of polygon are the same, if yes it compares the y coordinates in
        % order to decide where to insert the intersection point in polygon
        if Int(1,l) == polygon(1,i)
            if i ~=4 + insertions % for all but the last vertex of polygon
                % Tries to find a place between i-th and (i+1)-th vertex 
                if (Int(2,l) <= polygon(2,1+i) && Int(2,l) >= polygon(2,i)) ||...
                        (Int(2,l) >= polygon(2,1+i) && Int(2,l) <= polygon(2,i))
                    polygon = [polygon(:,1:i) Int(:,l) polygon(:,i+1:end)]; % insertion
                    flag1 = 1; % to brake the outer for loop and go to next intersection point
                    insertAfterPos = [insertAfterPos ,i]; % stores the location of insertion
                    insertions = insertions + 1; % counts the number of insertions
                    break;
                    
                end
            else % and now for the last vertex check if it fits between last and first vertex
                if Int(2,l) < polygon(2,4 + insertions) && Int(2,l) > polygon(2,1) ||...
                        Int(2,l) > polygon(2,4 + insertions) && Int(2,l) < polygon(2,1)
                    polygon = [polygon, Int(:,l)]; % insertion
                    flag1 = 1; % to brake the outer for loop and go to next intersection point
                    insertAfterPos = [insertAfterPos , 4 + insertions]; % stores the location of insertion
                    insertions = insertions + 1; % counts the number of insertions
                    break;
                    
                end
            end
            
            
            if flag1 == 1, flag1 = 0; break, end
         
        % Check if y coordinates of l-th intersection point and i-th vertex
        % of polygon are the same, if yes it compares the x coordinates in
        % order to decide where to insert the intersection point in polygon    
        elseif Int(2,l) == polygon(2,i)
            if i~= 4 + insertions % for all but the last vertex of polygon
                % Tries to find a place between i-th and (i+1)-th vertex 
                if (Int(1,l) <= polygon(1,1+i) && Int(1,l) >= polygon(1,i)) ||...
                        (Int(1,l) >= polygon(1,1+i) && Int(1,l) <= polygon(1,i))
                    polygon = [polygon(:,1:i) Int(:,l) polygon(:,i+1:end)]; % insertion
                    flag2 = 1; % to brake the outer for loop and go to next intersection point
                    insertAfterPos = [insertAfterPos ,i]; % stores the location of insertion
                    insertions = insertions + 1; % counts the number of insertions
                    break;

                end
            else % and now for the last vertex check if it fits between last and first vertex
                if Int(1,l) < polygon(1,4 + insertions) && Int(1,l) > polygon(1,1) ||...
                        Int(1,l) > polygon(1,4 + insertions) && Int(1,l) < polygon(1,1)
                    polygon = [polygon, Int(:,l)]; % insertion
                    flag2 = 1; % to brake the outer for loop and go to next intersection point
                    insertAfterPos = [insertAfterPos , 4 + insertions]; % stores the location of insertion
                    insertions = insertions + 1; % counts the number of insertions
                    break;

                end
            end
            if flag2 == 1, flag2 = 0; break, end
        end
    end
end

% The locations to cut the quad polygon when inserting the vertices of
% the control polygon
cut1 = insertAfterPos(:,1)+1;
cut2 = insertAfterPos(:,2)+1;

% Inserting the vertices of the control polygon/polygonised curve with 
%respect to the direction of rotation and generating the 2 subdomains 
% polygon 1 & 2. It is checked to which of the inserted intersection points
% the first coordiante of the array points is equal to and based on this,
% it is decided if the the aray points must be inserted in the current
% order or flipped.

if insertAfterPos(:,1) < insertAfterPos(:,2)
    
    if isequal(polygon(:,(cut2)),points(:,1))
        
        polygon1 = [polygon(:,(cut1):(cut2)), points(:,2:(length(points)-1))];
        polygon2 = [polygon(:,1:(cut1)), points_flip(:,2:(length(points)-1)), polygon(:,(cut2):end)];
   
    else
        
        polygon1 = [polygon(:,(cut1):(cut2)), points_flip(:,2:(length(points)-1))];
        polygon2 = [polygon(:,1:(cut1)), points(:,2:(length(points)-1)), polygon(:,(cut2):end)];

    end

else
    insertAfterPos(:,1) = insertAfterPos(:,1) + 1;
    cut1 = insertAfterPos(:,1)+1;
    if isequal(polygon(:,(cut1)),points(:,1))
        polygon1 = [polygon(:,(cut2):(cut1)), points(:,2:(length(points)-1))];
        polygon2 = [polygon(:,1:(cut2)), points_flip(:,2:(length(points)-1)), polygon(:,(cut1):end)];
    else
        polygon1 = [polygon(:,(cut2):(cut1)), points_flip(:,2:(length(points)-1))];
        polygon2 = [polygon(:,1:(cut2)), points(:,2:(length(points)-1)), polygon(:,(cut1):end)];
    end

end

end

