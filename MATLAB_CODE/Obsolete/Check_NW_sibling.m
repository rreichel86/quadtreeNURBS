function [extract_element]=Check_NW_sibling(Quadtree,i,l,Location,Loc_Current)
%  Check_NW_sibling function looks for the possible neighbours of a NW quad and
% takes out the coordinates of neighboring leafs.In NW leaf the neighboring
% element of East and South are of the same father as of NW.While the other
% two neighbors West and North are sons of the neighboring fathers.

% Input:
% Quadtree data
% Location: adressbook of all nodes
% Loc_Current: adress of given quad
% Output:
% extract_element:coordinates of neighboring edges that share boundary with
% current leaf

Father=Quadtree.Node{Quadtree.Parent(l(i)),1};% Pointer to father of actual quad
idx_Father=Quadtree.Parent(l(i));%it´s index

% Adress of south sibiling
if idx_Father==1;idx_S_Sib=Father{2}(2);else;idx_S_Sib=Father{11}(2);end
Loc_S_Sib=Location{idx_S_Sib};
% Checking the level of decomposition of it's south sibling
[nSons_S_Sib]= number_sonsext(Location,Loc_S_Sib);

% if South sibling have sons than the South quad siblings share the
% coordinates to current quad otherwise default -99
if nSons_S_Sib>1
    Father_South_Sib= Quadtree.Node{idx_S_Sib,1};
    idx_S_Sib_1=Father_South_Sib{11}(1);
    son_South_Sib=Quadtree.Node{idx_S_Sib_1,1};
    [Extracting_S_corxy]=son_South_Sib{10}(:,3);
else
    [Extracting_S_corxy]=zeros(2,1)-99;
end

% Adress of east sibiling
if idx_Father==1;idx_E_Sib=Father{2}(3);else;idx_E_Sib=Father{11}(3);end
Loc_E_Sib=Location{idx_E_Sib};
% Checking the level of decomposition of it's east sibling
[nSons_E_Sib]= number_sonsext(Location,Loc_E_Sib);

% if East sibling have sons than the East quad siblings share the
% coordinates to current quad otherwise default -99
if nSons_E_Sib>1
    Father_East_Sib= Quadtree.Node{idx_E_Sib,1};
    idx_E_Sib_1=Father_East_Sib{11}(1);
    son_East_Sib=Quadtree.Node{idx_E_Sib_1,1};
    [Extracting_E_corxy]=son_East_Sib{10}(:,1);
else
    [Extracting_E_corxy]=zeros(2,1)-99;
end

% Adress of eventual north neighbor quad,finding it's location
if length(Location{l(i)})>=1
    Loc_N_Sib=zeros(1,length(Location{l(i)}));
    copy=0;
    for j=length(Loc_N_Sib):-1:1
        if copy==1
            Loc_N_Sib(j)=Location{l(i)}(j);
        end
        if Location{l(i)}(j)==1 || Location{l(i)}(j)==3
            if copy==0;Loc_N_Sib(j)=Location{l(i)}(j)+1;end
        end
        if Location{l(i)}(j)==2 || Location{l(i)}(j)==4
            if copy==0;Loc_N_Sib(j)=Location{l(i)}(j)-1;end
            copy=1;
        end
    end
end

for j=1:length(Loc_N_Sib)
    if Loc_N_Sib(j)<1;Loc_N_Sib(1)=Loc_N_Sib(1)+4;
    end
    if Loc_N_Sib(j)>4;Loc_N_Sib(1)=Loc_N_Sib(1)-4;
    end
end
%if first index of current quad is 1(NW) and first index of
%North Neighbor quad is 2 then it means it's the top left quad
%so we don't need to have coordinates of north quad and it's
%same if current quad(NW) in 3(NE) and first index of
%North Neighbor quad is 4 than no need of north quad
if (Loc_Current(1)==1 && Loc_N_Sib(1)==2) || (Loc_Current(1)==3 && Loc_N_Sib(1)==4)
    [Extracting_N_corxy]=zeros(2,1)-99;
else
    % Finding if sibiling exists
    idx = cellfun('length',Location)==length(Loc_N_Sib);
    for j=1:length(idx)
        tf=isequal(Location{j},Loc_N_Sib);
        if tf==true;idx_N_Sib=j;
        end
    end
    
    
    [nSons_N_Sib]= number_sonsext(Location,Loc_N_Sib);
    % if North sibling have sons than the North quad sons share the
    % coordinates to current quad otherwise default -99
    if nSons_N_Sib>1
        Father_North_Sib= Quadtree.Node{idx_N_Sib,1};
        idx_N_Sib_1=Father_North_Sib{11}(4);
        son_North_Sib=Quadtree.Node{idx_N_Sib_1,1};
        [Extracting_N_corxy]=son_North_Sib{10}(:,1);
    else
        [Extracting_N_corxy]=zeros(2,1)-99;
        
    end
    
end

% Adress of eventual west neighbor quad,finding it's location
if length(Location{l(i)})>=1
    Loc_W_Sib=zeros(1,length(Location{l(i)}));
    copy=0;
    for j=length(Loc_W_Sib):-1:1
        if copy==1
            Loc_W_Sib(j)=Location{l(i)}(j);
        end
        if Location{l(i)}(j)==1 || Location{l(i)}(j)==2
            if copy==0;Loc_W_Sib(j)=Location{l(i)}(j)+2;end
        end
        if Location{l(i)}(j)==3 || Location{l(i)}(j)==4
            if copy==0;Loc_W_Sib(j)=Location{l(i)}(j)-2;end
            copy=1;
        end
    end
end

for j=1:length(Loc_W_Sib)
    if Loc_W_Sib(j)<1;Loc_W_Sib(1)=Loc_W_Sib(1)+4;
    end
    if Loc_W_Sib(j)>4;Loc_W_Sib(1)=Loc_W_Sib(1)-4;
    end
end
%if first index of current quad is 1(NW) and first index of
%West Neighbor quad is 3 then it means it's the top left quad
%so we don't need to have coordinates of west quad and it's
%same if current quad(NW) in 2(SW) and first index of
%West Neighbor quad is 4 than no need of west quad
if (Loc_Current(1)==1 && Loc_W_Sib(1)==3) || (Loc_Current(1)==2 && Loc_W_Sib(1)==4)
    [Extracting_W_corxy]=zeros(2,1)-99;
else
    % Finding index of west neighbor
    idx = cellfun('length',Location)==length(Loc_W_Sib);
    for j=1:length(idx);tf=isequal(Location{j},Loc_W_Sib);
        if tf==true;idx_W_Sib=j;
        end
    end
    
    [nSons_W_Sib]= number_sonsext(Location,Loc_W_Sib);
    if nSons_W_Sib>1
        % if West Neighbor have sons than the West quad sons share the
        % coordinates to current quad otherwise default -99
        Father_West_Sib= Quadtree.Node{idx_W_Sib,1};
        idx_W_Sib_1=Father_West_Sib{11}(4);
        son_West_Sib=Quadtree.Node{idx_W_Sib_1,1};
        [Extracting_W_corxy]=son_West_Sib{10}(:,3);
    else
        [Extracting_W_corxy]=zeros(2,1)-99;
    end
    
end
extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];

end



