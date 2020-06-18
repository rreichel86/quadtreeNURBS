 function [Extracting_S_corxy,Extracting_E_corxy,Extracting_N_corxy,Extracting_W_corxy]=Check_NW_sibling(Quadtree,i,l,...
        Location,Loc_Current);  
        
        Father=Quadtree.Node{Quadtree.Parent(l(i)),1};
        idx_Father=Quadtree.Parent(l(i));%it´s index
       
      % Adress of south sibiling
        if idx_Father==1;idx_S_Sib=Father{2}(2);else;idx_S_Sib=Father{11}(2);end
        Loc_S_Sib=Location{idx_S_Sib};
        % Checking it's level of decomposition
        [nSons_S_Sib]= number_sonsext(Location,Loc_S_Sib);
        
        % if level of decompisition differ in more than one, decompose
        if nSons_S_Sib>1
            % Decompose_balance functions prepares input for decompose function 
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
        % Checking it's level of decomposition
        [nSons_E_Sib]= number_sonsext(Location,Loc_E_Sib);
        
        % if level of decompisition differ in more than one, decompose
        if nSons_E_Sib>1
            Father_East_Sib= Quadtree.Node{idx_E_Sib,1};
            idx_E_Sib_1=Father_East_Sib{11}(1);
            son_East_Sib=Quadtree.Node{idx_E_Sib_1,1};
            [Extracting_E_corxy]=son_East_Sib{10}(:,1);
             else
            [Extracting_E_corxy]=zeros(2,1)-99;
        end
        
        % Adress of eventual north sibiling, findig if it exists and level
        % of decompoition
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
            
            if (Loc_Current(1)==1 & Loc_N_Sib(1)==2) | (Loc_Current(1)==3 & Loc_N_Sib(1)==4)
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
            % if level of decompisition differ in more than one, decompose
            if nSons_N_Sib>1
            % Decompose_balance functions prepares input for decompose function
                Father_North_Sib= Quadtree.Node{idx_N_Sib,1};
                idx_N_Sib_1=Father_North_Sib{11}(4);
                son_North_Sib=Quadtree.Node{idx_N_Sib_1,1};
                [Extracting_N_corxy]=son_North_Sib{10}(:,1);
            else
                [Extracting_N_corxy]=zeros(2,1)-99;
                
        end
          
            end
        
        % Adress of eventual west sibiling, findig if it exists and level
        % of decompoition
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
            
          if (Loc_Current(1)==1 & Loc_W_Sib(1)==3) | (Loc_Current(1)==2 & Loc_W_Sib(1)==4)
            [Extracting_W_corxy]=zeros(2,1)-99;
            else      
            % Finding if sibiling exists
            idx = cellfun('length',Location)==length(Loc_W_Sib);
            for j=1:length(idx);tf=isequal(Location{j},Loc_W_Sib);
                if tf==true;idx_W_Sib=j;
                end;
            end
    
                [nSons_W_Sib]= number_sonsext(Location,Loc_W_Sib);
%                 % if level of decompisition differ in more than one, decompose
                if nSons_W_Sib>1
%                     % Decompose_balance functions prepares input for decompose function   
               Father_West_Sib= Quadtree.Node{idx_W_Sib,1};
                idx_W_Sib_1=Father_West_Sib{11}(4);
                son_West_Sib=Quadtree.Node{idx_W_Sib_1,1};
                [Extracting_W_corxy]=son_West_Sib{10}(:,3);
            else
                [Extracting_W_corxy]=zeros(2,1)-99;
                end
                
          end
              
 end



