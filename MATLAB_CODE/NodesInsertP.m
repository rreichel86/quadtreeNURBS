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
%                               nodes, where the first three entries
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
% seedingPoints = [isec,isec0,idxLeaf,xcoor,ycoor]
%
%                              isec  - new section number of the (un)qualified section
%                              isec0 - old section number lof the (un)qualified section 
%                              idxLeaf - number of Leaf
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
pgrad_splitt = max(seedingPoints_splitt(:,6));
pgrad_merge = max(seedingPoints_merge(:,6));

if pgrad_splitt + 1 <= pgrad_merge 
    sections = [sections,zeros(length(sections(:,1)),pgrad_merge-2)];
else
    sections = [sections,zeros(length(sections(:,1)),pgrad_splitt-1)];
end

%% insert nodes in P-direction
num_seedingPoints = length(seedingPoints_splitt(:,1));
sections = [sections,zeros(length(sections(:,1)),ep-1)];
ninode = ep - 1;
for isp = 1: num_seedingPoints           
    isec0 = seedingPoints_splitt(isp,2);
    ikv = sections(isec0,4); 
    coor_nsec = []; %coor_matrix contains the information of all nodes including scalling center for each chosed section

    % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
    if ikv == 0
        if sections(isec0,6) == 3 %check if this section has been treated in secN

            nodesNum = sections(isec0,7:9); %[node_1,node_2,scaling center]
            for ii = 1:length(nodesNum)
                number = nodesNum(1,ii);
                coor_nsec0 = coor(number,:);%extract nodes coordinates of each selected section
                coor_nsec = [coor_nsec;coor_nsec0];
            end
            %the coordinate of the inserted node in p-direction           
            coor_nodes_p = getLocation(coor_nsec(1,2),coor_nsec(2,2),coor_nsec(1,3),coor_nsec(2,3),ep);
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_nodes_p];
            % update the section matrix 
            sections(isec0,6) = 3 + ninode;
            sections(isec0,end-1) = sections(isec0,9);%shift of the scaling center
            sections(isec0,end-2) = sections(isec0,8);%shift of the structral points
            sections(isec0,8:7+ninode) = (nnode+1:nnode+ninode);

            %update the ord matrix       
            ord(isec0,2) = ep;

            idx_isecN = find(secN_splitt(:,1) == isec0);
            isecN = secN_splitt(idx_isecN,2);
            if isecN ~= 0 
                sections(isecN,6) = 3 + ninode; 
                sections(isecN,end-1) = sections(isecN,9);%shift of the scaling center
                sections(isecN,end-2) = sections(isecN,8);%shift of the structural point
                sections(isecN,8:7+ninode) = rot90(sections(isec0,8:7+ninode),2);
                ord(isecN,2) = ep;
            end
            nnode = nnode + ninode;
        end        
    end
end
            

end    