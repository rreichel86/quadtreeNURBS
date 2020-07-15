function [extract_element]=Check_SW_sibling(Quadtree,i,l,Location,Loc_Current)
%  Check_SW_sibling function looks for the possible neighbours of a SW quad and
% takes out the coordinates of neighboring leafs.In SW leaf the neighboring
% element of East and North are of the same father as of SW.While the other
% two neighbors West and South are sons of the neighboring fathers.

% Input:
% Quadtree data
% Location: adressbook of all nodes
% Loc_Current: adress of given quad
% Output:
% extract_element:coordinates of neighboring edges that share boundary with
% current leaf

% Adress of eventual south neighbor quad,finding it's location
if length(Location{l(i)})>=1
    Loc_S_Sib=zeros(1,length(Location{l(i)}));
    copy=0;
    for j=length(Loc_S_Sib):-1:1
        if copy==1
            Loc_S_Sib(j)=Location{l(i)}(j);
        end
        if Location{l(i)}(j)==2 || Location{l(i)}(j)==4
            if copy==0;Loc_S_Sib(j)=Location{l(i)}(j)-1;end
        end
        if Location{l(i)}(j)==1 || Location{l(i)}(j)==3
            if copy==0;Loc_S_Sib(j)=Location{l(i)}(j)+1;end
            copy=1;
        end
    end
end

for j=1:length(Loc_S_Sib)
    if Loc_S_Sib(j)<1
        Loc_S_Sib(1)=Loc_S_Sib(1)+4;
    end
    if Loc_S_Sib(j)>4
        Loc_S_Sib(1)=Loc_S_Sib(1)-4;
    end
end
%if first index of current quad is 2(SW) and first index of
%south Neighbor quad is 1 then it means it's the bottom left quad
%so we don't need to have coordinates of south quad and it's
%same if current quad(SW) in 4(SW) and first index of
%South Neighbor quad is 3 than no need of south quad
if (Loc_Current(1)==2 && Loc_S_Sib(1)==1) || (Loc_Current(1)==4 && Loc_S_Sib(1)==3)
    [Extracting_S_corxy]=zeros(2,1)-99;
else
    % Finding if sibiling exists
    idx = cellfun('length',Location)==length(Loc_S_Sib);
    for j=1:length(idx)
        tf=isequal(Location{j},Loc_S_Sib);
        if tf==true;idx_S_Sib=j;
        end
    end
    
    [nSons_S_Sib]= number_sonsext(Location,Loc_S_Sib);
    % if South sibling have sons than the South quad sons share the
    % coordinates to current quad otherwise default -99
    if nSons_S_Sib>1
        Father_South_Sib= Quadtree.Node{idx_S_Sib,1};
        idx_S_Sib_1=Father_South_Sib{11}(1);
        son_South_Sib=Quadtree.Node{idx_S_Sib_1,1};
        [Extracting_S_corxy]=son_South_Sib{10}(:,3);
    else
        [Extracting_S_corxy]=zeros(2,1)-99;
    end
end



Father=Quadtree.Node{Quadtree.Parent(l(i)),1};% Pointer to father of actual quad
idx_Father=Quadtree.Parent(l(i));%it´s index


% Adress of east sibiling
if idx_Father==1;idx_E_Sib=Father{2}(4);else;idx_E_Sib=Father{11}(4);end
Loc_E_Sib=Location{idx_E_Sib};
% Checking the level of decomposition of it's east sibling
[nSons_E_Sib]= number_sonsext(Location,Loc_E_Sib);

% if East sibling have sons than the East quad siblings share the
% coordinates to current quad otherwise default -99
if nSons_E_Sib>1
    Father_East_Sib= Quadtree.Node{idx_E_Sib,1};
    idx_E_Sib_1=Father_East_Sib{11}(2);
    son_East_Sib=Quadtree.Node{idx_E_Sib_1,1};
    [Extracting_E_corxy]=son_East_Sib{10}(:,4);
else
    [Extracting_E_corxy]=zeros(2,1)-99;
end


% Adress of north sibiling
if idx_Father==1;idx_N_Sib=Father{2}(1);else;idx_N_Sib=Father{11}(1);end
Loc_N_Sib=Location{idx_N_Sib};
% Checking the level of decomposition of it's north sibling
[nSons_N_Sib]= number_sonsext(Location,Loc_N_Sib);

% if North sibling have sons than the North quad siblings share the
% coordinates to current quad otherwise default -99
if nSons_N_Sib>1
    
    Father_North_Sib= Quadtree.Node{idx_N_Sib,1};
    idx_N_Sib_1=Father_North_Sib{11}(2);
    son_North_Sib=Quadtree.Node{idx_N_Sib_1,1};
    [Extracting_N_corxy]=son_North_Sib{10}(:,2);
else
    [Extracting_N_corxy]=zeros(2,1)-99;
end



% Adress of eventual west sibiling, findig it's location
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
    if Loc_W_Sib(j)<1
        Loc_W_Sib(1)=Loc_W_Sib(1)+4;
    end
    if Loc_W_Sib(j)>4;Loc_W_Sib(1)=Loc_W_Sib(1)-4;
    end
end
%if first index of current quad is 1(SW) and first index of
%West Neighbor quad is 3 then  we don't need to have
%coordinates of west quad and it's same if current quad(SW) in
%2(SW) and first index of West Neighbor quad is 4 than no need
%of West quad
if (Loc_Current(1)==1 && Loc_W_Sib(1)==3) || (Loc_Current(1)==2 && Loc_W_Sib(1)==4)
    [Extracting_W_corxy]=zeros(2,1)-99;
else
    % Finding if sibiling exists
    idx = cellfun('length',Location)==length(Loc_W_Sib);%finding if sibiling exists
    for j=1:length(idx)
        tf=isequal(Location{j},Loc_W_Sib);
        if tf==true;idx_W_Sib=j;
        end
    end
    
    
    [nSons_W_Sib]= number_sonsext(Location,Loc_W_Sib);
    % if West sibling have sons than the West quad sons share the
    % coordinates to current quad otherwise default -99
    if nSons_W_Sib>1
        Father_West_Sib= Quadtree.Node{idx_W_Sib,1};
        idx_W_Sib_1=Father_West_Sib{11}(3);
        son_West_Sib=Quadtree.Node{idx_W_Sib_1,1};
        [Extracting_W_corxy]=son_West_Sib{10}(:,2);
    else
        [Extracting_W_corxy]=zeros(2,1)-99;
    end
end
extract_element=[Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy];
end






