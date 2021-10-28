function [coor,nnode,sections,ord] = NodesInsertP(ep,nnode,coor,sections,seedingPoints_splitt,ord,secN) 
% Insert nodes in p_direction

%INPUT:
%ep = elevated pgrad 
%eq = elevated qgrad 
%
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                                    3 - inserted nodes
%                              which_region: region number
%                              inside_region: 0 - at the boundary (NURBS
%                                                         curve boundary?)
%                                             1 - inside (inside of hole?)
%                                            -1 - outside (out of hole?)
%
% sections -------------------- sections connectivity matrix as nsec-tupel of 
%                               nodes, where the first three entries
%                               isec - section number
%                               ikv - knot vector number
%                               iel - element number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, idxLeaf, ikv, iel,region, nsec, node_1,...,node_nsec]
%
%
% seedingPoints_splitt = [isec,isec0,idxLeaf,xcoor,ycoor,c]
%                         
%                       isec  - new section number of the unqualified section
%                       isec0 - old section number lof the unqualified section 
%                       idxLeaf - number of Leaf
%                       x/y_coor - x/y coordinate of the scaling center
%                                  of the unqualified sections
%                       c  -   error_measure
%
% ord = [isec,pgrad,qgrad]
%
% secN = [isec, isecN]
%
% polyElmts -------------------- relate sections and polygonal elements
% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
%
%
%
%OUTPUT: coor, nnode, sections, ord                   


num_seedingPoints = length(seedingPoints_splitt(:,1));

%% insert nodes in P-direction
sections = [sections,zeros(length(sections(:,1)),ep-1)];
ninode = ep - 1;
for isp = 1: num_seedingPoints           
    isec0 = seedingPoints_splitt(isp,2);
    ikv = sections(isec0,3); 
    coor_nsec = []; %coor_matrix contains the information of all nodes including scalling center for each chosed section
    % sections = [isec, idxLeaf, ikv, iel,region, nsec, node_1,...,node_nsec]
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

            idx_isecN = find(secN(:,1) == isec0);
            isecN = secN(idx_isecN,2);
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