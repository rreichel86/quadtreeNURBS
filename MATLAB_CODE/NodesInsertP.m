function [coor,nnode,sections,ord] = NodesInsertP(nnode,coor,sections,ord,seedingPoints_splitt,secN_splitt,seedingPoints_merge,secN_merge) 
% Insert nodes in p_direction

%INPUT:
%
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 
%
% sections -------------------- sections connectivity matrix as nsec-tupel of 
%                               nodes, where the first six entries
%                               isec - section number
%                               ipoly - polygonal element number 
%                               idxLeaf - index of leaf 
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
%
%
% ord = [isec,pgrad,qgrad]
%
% seedingPoints = [isec,isec0,idxLeaf,xcoor,ycoor,pgrad,qgrad]
%
%                              isec  - new section number of the (un)qualified section
%                              isec0 - old section number lof the (un)qualified section 
%                              idxLeaf - number of Leaf
%                              knov  - knot vector
%                              x/y_coor - x/y coordinate of the scaling center
%                                         of the unqualified sections
%                              p-/qgrad  - p-/qgrad from last calculation
%
% secN_splitt = [isec_splitt, isecN_splitt] ---  section number of unqualified section and its neighbour section 
% secN_merge  = [isec_merge, isecN_merge]  ---   section number of qualified section and its neighbour section
%
%
%
%OUTPUT: coor, nnode, sections, ord       

% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                                    3 - inserted node
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 


%% extend sections matrix

% max. pgrad of all unqualified sections from last calculation
pgrad_splitt = max(seedingPoints_splitt(:,7));

if isempty(seedingPoints_merge) == 0
    
    % max. pgrad of all qualified sections from last calculation
    pgrad_merge = max(seedingPoints_merge(:,7));
else
    pgrad_merge = 1;
end 

% number of sections 
numsec = length(sections(:,1));

if pgrad_splitt + 1 <= pgrad_merge 
    sections = [sections,zeros(numsec,pgrad_merge-2)];
else
    sections = [sections,zeros(numsec,pgrad_splitt-1)];
end

%% insert nodes in P-direction for unqualified sections

numSeedingPoints_splitt = length(seedingPoints_splitt(:,1));

for isp = 1: numSeedingPoints_splitt           
    isec0 = seedingPoints_splitt(isp,2); %(old)section number
    pgrad = seedingPoints_splitt(isp,7); %pgrad from last calculation 
    ep = pgrad + 1; %elevated pgrad
    ninode = ep - 1; %number of inserted nodes
    ikv = sections(isec0,4); %knot vector
    
    if ninode == 0
        break
    end
    
    coor_StrucNodes = []; %coor_matrix contains the information of all nodes including scalling center for each chosed section
    % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
    if ikv == 0
        if sections(isec0,6) == 3 %check if this section has been treated in secN

            StrucNodes = sections(isec0,7:9); %[node_1,node_2,scaling center]
            for ii = 1:length(StrucNodes)
                number = StrucNodes(1,ii);
                coor_StrucNodes0 = coor(number,:);%extract nodes coordinates of each selected section
                coor_StrucNodes = [coor_StrucNodes;coor_StrucNodes0];
            end
            %the coordinate of the inserted node in p-direction           
            coor_nodes_p = getLocation(coor_StrucNodes(1,2),coor_StrucNodes(2,2),coor_StrucNodes(1,3),coor_StrucNodes(2,3),ep);
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_nodes_p];
            % update the section matrix 
            nsec_new = 3 + ninode;
            sections(isec0,6) = nsec_new;
            sections(isec0,6+nsec_new) = sections(isec0,9);%shift of the scaling center
            sections(isec0,5+nsec_new) = sections(isec0,8);%shift of the structral points
            sections(isec0,8:7+ninode) = (nnode+1:nnode+ninode);

            %update the ord matrix       
            ord(isec0,2) = ep;

            idx_isecN = find(secN_splitt(:,1) == isec0);
            isecN = secN_splitt(idx_isecN,2);
            if isecN ~= 0
                nsecN_new = 3 + ninode;
                sections(isecN,6) = nsecN_new;                 
                sections(isecN,6+nsecN_new) = sections(isecN,9);%shift of the scaling center
                sections(isecN,5+nsecN_new) = sections(isecN,8);%shift of the structural point
                sections(isecN,8:7+ninode) = rot90(sections(isec0,8:7+ninode),2);
                ord(isecN,2) = ep;
            end
            nnode = nnode + ninode;
        end        
    end
end


%% insert nodes in P-direction for qualified sections

if isempty(seedingPoints_merge) == 1
    return
end

numSeedingPoints_merge = length(seedingPoints_merge(:,1));

for isp = 1: numSeedingPoints_merge           
    isec0 = seedingPoints_merge(isp,2); %(old)section number
    pgrad = seedingPoints_merge(isp,7); %pgrad from last calculation
    ninode = pgrad - 1; %number of inserted nodes
    ikv = sections(isec0,4); %knot vector
    
    if ninode == 0
        break
    end

    coor_StrucNodes = []; %coordinate of structural nodes including scalling center for each chosed section
    % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
    
    %check if current section is the neighbor section of one unqualified section
    if ismember(isec0,secN_splitt(:,2)) ~= 1 %no        
        if ikv == 0
            if sections(isec0,6) == 3 %check if this section has been treated in secN_merge

                StrucNodes = sections(isec0,7:9); %[node_1,node_2,scaling center]
                for ii = 1:length(StrucNodes)
                    number = StrucNodes(1,ii);
                    coor_StrucNodes0 = coor(number,:);%extract nodes coordinates of each selected section
                    coor_StrucNodes = [coor_StrucNodes;coor_StrucNodes0];
                end
                %the coordinate of the inserted node in p-direction           
                coor_nodes_p = getLocation(coor_StrucNodes(1,2),coor_StrucNodes(2,2),coor_StrucNodes(1,3),coor_StrucNodes(2,3),pgrad);
                coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                coor(nnode+1:nnode+ninode,4) = 1;  %weight
                coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
                coor(nnode+1:nnode+ninode,2:3) =[coor_nodes_p];
                % update the section matrix
                nsec_new = 3 + ninode;
                sections(isec0,6) = nsec_new;
                sections(isec0,6+nsec_new) = sections(isec0,9);%shift of the scaling center
                sections(isec0,5+nsec_new) = sections(isec0,8);%shift of the structral points
                sections(isec0,8:7+ninode) = (nnode+1:nnode+ninode);

                %update the ord matrix       
                ord(isec0,2) = pgrad;

                idx_isecN = find(secN_merge(:,1) == isec0);
                isecN = secN_merge(idx_isecN,2);
                nsecN = sections(idx_isecN,6);
                if isecN ~= 0 

                    nsecN_new = 3 + ninode;
                    sections(isecN,6) = nsecN_new; 
                    sections(isecN,6+nsecN_new) = sections(isecN,9);%shift of the scaling center
                    sections(isecN,5+nsecN_new) = sections(isecN,8);%shift of the structural point
                    sections(isecN,8:7+ninode) = rot90(sections(isec0,8:7+ninode),2);
                    ord(isecN,2) = pgrad;
                end
                nnode = nnode + ninode;
            end        
        end
    end
end


end    