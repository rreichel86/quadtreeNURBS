function [extract_element] = extracting_element(Quadtree,l,i);
%Function to get the coordinates of neighboring element that are
%sharing the boundary
% Input:
% Quadtree data
%Quadleaf indices in QuadtreeNode data
% Definition of the NURBS
% Output:
% extract_element

% The function Location_Quads creates an adressbook of all nodes in 1:4
% format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
Location=Location_Quads(Quadtree);


    % The last two elements of its direction determine if 
    % the quad is NW,SW,NE,SE
    P=[Quadtree.Node{l(i),1}{2,1}(end-1) Quadtree.Node{l(i),1}{2,1}(end)];
    Loc_Current=Location{l(i)};%its location in 1:4 format
   

    % if we have a NW Quad. Extract the coordinates of neighbouring
    % leafs if they exist otherwise it would be -99
    if P==[1 1]
    [extract_element]=Check_NW_sibling(Quadtree,i,l,Location,Loc_Current);   
    end    
    
    %if we have a SW Quad.Extract the coordinates of neighbouring
    % leafs if they exist otherwise it would be -99
    if P==[2 1]
    [extract_element]=Check_SW_sibling(Quadtree,i,l,Location,Loc_Current);
    end
    
    %if we have a NE Quad.Extract the coordinates of neighbouring
    % leafs if they exist otherwise it would be -99
    if P==[1 2]
    [extract_element]=Check_NE_sibling(Quadtree,i,l,Location,Loc_Current);
    end
    
    %if we have a SE Quad.Extract the coordinates of neighbouring
    % leafs if they exist otherwise it would be -99
    if P==[2 2]
    [extract_element]=Check_SE_sibling(Quadtree,i,l,Location,Loc_Current);
    end 
    
end   