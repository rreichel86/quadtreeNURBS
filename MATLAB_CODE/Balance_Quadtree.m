function [Quadtree] = Balance_Quadtree(Quadtree,NURBS,controlPoints,...
    knots, weights, degree,Boundary)
% Balance Quadtree function obtains a balanced Quadtree out of an
% unbalanced tree. It goes through the nodes of the tree, checking how
% manys subdivions it has, and comparing it with the neighbours. If any
% quad is unbalanced, it calls the decomposition function. This is done
% recursively done until the tree is balanced  
%
% Input:
% Quadtree data
% Definition of the NURBS
% Output:
% Balanced Quadtree

% Variables initialization
i=2;
update_Location=0;

% The function Location_Quads creates an adressbook of all nodes in 1:4
% format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
Location=Location_Quads(Quadtree);


while i<=length(Quadtree.Node)
    % loop over the nodes of the tree. The last two elements of its 
    % direction determine if the quad is NW,SW,NE,SE
    P=[Quadtree.Node{i,1}{2,1}(end-1) Quadtree.Node{i,1}{2,1}(end)];
    Loc_Current=Location{i};%its location in 1:4 format
    
    % Obtaining ammount of sons' generations (level of decomposition) of 
    % current quad. This will be compared with it's neighbours. The
    % function number_sons takes the current location and the adressbook
    % 'Location' and gives as output the number of sons' generations
    [nSons]= number_sons(Location,Loc_Current);


    % if we have a NW Quad. Check_ functions compare the number of sons' 
    % generations for a current NW quad with it's neighbours and calls 
    % eventually the decomposition function
    if P==[1 1]
    [Quadtree,update_Location]=Check_NW(Quadtree, NURBS, controlPoints,...
      knots, weights, degree,i,Location,nSons,Loc_Current,Boundary);   
    end    
    %if we have a SW Quad
    if P==[2 1]
    [Quadtree,update_Location]=Check_SW(Quadtree, NURBS, controlPoints,...
      knots, weights, degree,i,Location,nSons,Loc_Current,Boundary);
    end
    %if we have a NE Quad
    if P==[1 2]
    [Quadtree,update_Location]=Check_NE(Quadtree, NURBS, controlPoints,...
      knots, weights, degree,i,Location,nSons,Loc_Current,Boundary);
    end
    %if we have a SE Quad
    if P==[2 2]
    [Quadtree,update_Location]=Check_SE(Quadtree, NURBS, controlPoints,...
      knots, weights, degree,i,Location,nSons,Loc_Current,Boundary);
    end 
    
if update_Location==1
    Location=Location_Quads(Quadtree);
    i=1;
end
i=i+1;    
end
end