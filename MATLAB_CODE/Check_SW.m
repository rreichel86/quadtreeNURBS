function [Quadtree,update_Location]=Check_SW(Quadtree, NURBS, controlPoints,...
    knots, weights, degree,i,Location,nSons,Loc_Current,Boundary)
% Check_SW function looks for the possible neighbours of a SW quad and
% compares their level of decomponsitions. If they differ in more than one
% level, the 2:1 rule is not fulfilled, therefore it calls decomposition
%
% Input:
% Quadtree data
% Definition of the NURBS
% Location: adressbook of all nodes
% nSons: ammount of sons' generations of given quad
% Loc_Current: adress of given quad
% Output:
% Quadtree after eventual decomposition
% update_Location: flag for updating adressbook

% Variables initialization
update_Location=0;
check=0;

% Pointer to eventual SW child
Location_Son=[Location{i} 1];
idx = cellfun('length',Location)==length(Location_Son);
for j=1:length(idx)
    tf=isequal(Location{j},Location_Son);if tf==true;check=1;end
end

% If the current quad already has children (check==1) we can't decompose it
if check == 0
    for iii=1:1
        
        % Pointer to father of actual quad
        Father=Quadtree.Node{Quadtree.Parent(i),1};
        idx_Father=Quadtree.Parent(i);
        check=0;
        
        % Adress of east sibiling
        if idx_Father==1;idx_E_Sib=Father{2}(4);else;idx_E_Sib=Father{11}(4);end
        Loc_E_Sib=Location{idx_E_Sib};
        % Checking it's level of decomposition
        [nSons_E_Sib]= number_sons(Location,Loc_E_Sib);
        
        % if level of decompisition differ in more than one, decompose
        if nSons_E_Sib>nSons+1
            % Decompose_balance functions prepares input for decompose function
            [Quadtree]=Decompose_balance(Quadtree,NURBS,controlPoints, ...
                knots, weights, degree,i,Loc_Current,idx_Father,Boundary);
            update_Location=1;
        end
        
        if update_Location==1;break;end
        check=0;
        
        % Adress of eventual south sibiling, findig if it exists and level
        % of decompoition
        if length(Location{i})>1
            Loc_S_Sib=zeros(1,length(Location{i}));
            copy=0;
            for j=length(Loc_S_Sib):-1:1
                if copy==1
                    Loc_S_Sib(j)=Location{i}(j);
                end
                if Location{i}(j)==2 || Location{i}(j)==4
                    if copy==0;Loc_S_Sib(j)=Location{i}(j)-1;end
                end
                if Location{i}(j)==1 || Location{i}(j)==3
                    if copy==0;Loc_S_Sib(j)=Location{i}(j)+1;end
                    copy=1;
                end
            end
            copy=0;
            
            for j=1:length(Loc_S_Sib)
                if Loc_S_Sib(j)<1
                    Loc_S_Sib(1)=Loc_S_Sib(1)+4;
                end
                if Loc_S_Sib(j)>4
                    Loc_S_Sib(1)=Loc_S_Sib(1)-4;
                end
            end
            
            % Finding if sibiling exists
            idx = cellfun('length',Location)==length(Loc_S_Sib);
            for j=1:length(idx)
                tf=isequal(Location{j},Loc_S_Sib);
                if tf==true;idx_S_Sib=j;check=1;
                end
            end
            
            if check==1
                [nSons_S_Sib]= number_sons(Location,Loc_S_Sib);
                % if level of decompisition differ in more than one, decompose
                if nSons_S_Sib>nSons+1
                    [Quadtree]=Decompose_balance(Quadtree,NURBS,controlPoints, ...
                        knots, weights, degree,i,Loc_Current,idx_Father,Boundary);
                    update_Location=1;
                end
            end
        end
        
        if update_Location==1;break;end
        check=0;
        
        
        % Adress of eventual west sibiling, findig if it exists and level
        % of decompoition
        if length(Location{i})>1
            Loc_W_Sib=zeros(1,length(Location{i}));
            copy=0;
            for j=length(Loc_W_Sib):-1:1
                if copy==1
                    Loc_W_Sib(j)=Location{i}(j);
                end
                if Location{i}(j)==1 || Location{i}(j)==2
                    if copy==0;Loc_W_Sib(j)=Location{i}(j)+2;end
                end
                if Location{i}(j)==3 || Location{i}(j)==4
                    if copy==0;Loc_W_Sib(j)=Location{i}(j)-2;end
                    copy=1;
                end
            end
            copy=0;
            
            for j=1:length(Loc_W_Sib)
                if Loc_W_Sib(j)<1
                    Loc_W_Sib(1)=Loc_W_Sib(1)+4;
                end
                if Loc_W_Sib(j)>4;Loc_W_Sib(1)=Loc_W_Sib(1)-4;
                end
            end
            
            % Finding if sibiling exists
            idx = cellfun('length',Location)==length(Loc_W_Sib);%finding if sibiling exists
            for j=1:length(idx)
                tf=isequal(Location{j},Loc_W_Sib);
                if tf==true;idx_W_Sib=j;check=1;
                end
            end
            
            if check==1
                [nSons_W_Sib]= number_sons(Location,Loc_W_Sib);
                % if level of decompisition differ in more than one, decompose
                if nSons_W_Sib>nSons+1
                    % Decompose_balance functions prepares input for decompose function
                    [Quadtree]=Decompose_balance(Quadtree,NURBS,controlPoints, ...
                        knots, weights, degree,i,Loc_Current,idx_Father,Boundary);
                    update_Location=1;
                end
            end
        end
        
        if update_Location==1;break;end
        check=0;
        
        
        % Adress of north sibiling
        if idx_Father==1;idx_N_Sib=Father{2}(1);else;idx_N_Sib=Father{11}(1);end
        Loc_N_Sib=Location{idx_N_Sib};
        % Checking it's level of decomposition
        [nSons_N_Sib]= number_sons(Location,Loc_N_Sib);
        
        % if level of decompisition differ in more than one, decompose
        if nSons_N_Sib>nSons+1
            % Decompose_balance functions prepares input for decompose function
            [Quadtree]=Decompose_balance(Quadtree,NURBS,controlPoints, ...
                knots, weights, degree,i,Loc_Current,idx_Father,Boundary);
            update_Location=1;
        end
    end
end
end



