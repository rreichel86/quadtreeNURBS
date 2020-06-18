function [extract_element] = extracting_element(Quadtree,l,i);


% The function Location_Quads creates an adressbook of all nodes in 1:4
% format (11 = 1, 21 = 2, 12 = 3, 22 = 4)
Location=Location_Quads(Quadtree);


    % loop over the nodes of the tree. The last two elements of its 
    % direction determine if the quad is NW,SW,NE,SE
    P=[Quadtree.Node{l(i),1}{2,1}(end-1) Quadtree.Node{l(i),1}{2,1}(end)];
    Loc_Current=Location{l(i)};%its location in 1:4 format
   

    % if we have a NW Quad. Check_ functions compare the number of sons' 
    % generations for a current NW quad with it's neighbours and calls 
    % eventually the decomposition function
    if P==[1 1]
    [Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy]=Check_NW_sibling(Quadtree,i,l,...
        Location,Loc_Current);   
  extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];
    end    
    %if we have a SW Quad
    if P==[2 1]
    [Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy]=Check_SW_sibling(Quadtree,i,l,...
        Location,Loc_Current);
  extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];
    end
    %if we have a NE Quad
    if P==[1 2]
    [Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy]=Check_NE_sibling(Quadtree,i,l,...
        Location,Loc_Current);
  extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];

    end
    %if we have a SE Quad
    if P==[2 2]
    [Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy]=Check_SE_sibling(Quadtree,i,l,...
        Location,Loc_Current);
  extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];

    end 
    
end   