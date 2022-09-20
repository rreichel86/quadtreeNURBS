function [coor,maxnsec,nnode,sections,ord]=elevateOrder_hRef_epeq(Quadtree,nnode,coor,sections,ord,polyElmts,connectivity,ep,eq)
% elevate order in p-direction and/or q-direction
%
%INPUT:

% Quadtree
% nnode = number of nodes 
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
%                               idxLeaf - number of leaf
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
%
%
% ord = [isec,pgrad,qgrad]
%
% polyElmts -------------------- relate sections and polygonal elements
% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
%
% connectivity --------------- elements connectivity matrix as nel-tupel of 
%                              nodes, where the first three entries
%                              iel - element number
%                              ikv - knot vector number
%                              idxLeaf - index of Leaf
%                              which_region - region number
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]

% seedingPoints = [isec,isec0,idxLeaf,xcoor,ycoor,pgrad,qgrad]
%                         
%                       isec  --- new section number of the (un)qualified section
%                       isec0 --- old section number lof the (un)qualified section 
%                       idxLeaf -- index of Leaf
%                       x/y_coor - x/y coordinate of the scaling center
%                       p-/qgrad  -- p-/qgrad from last calculation
%
%
%
%OUTPUT: coor, maxnsec, nnode, sections, ord 
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                                    3 - inserted nodes
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 


%% get the neighbours of sections
numSec = size(sections,1);

secNQ = cell(numSec,1); %array for all section numbers of neighbour quad
%secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

secN = []; %array for section numbers of the neighbour section
%secN = [idx_sec,isecN]


for isec = 1: numSec
    ikvo = sections(isec,4);    
    idxLeaf_sec = sections(isec, 3); %idxLeaf of sections
      
    refLeaf_sec = Quadtree.Node{idxLeaf_sec,1}{2,1}(1:end);
    secNQ{isec}(1,1) = isec;
   
    % search for current Quad neighbours
    % Loop over directions:
    % 1 - West
    % 2 - South
    % 3 - East
    % 4 - North
    idxNQa = [];
    for dir = 1:4 % determine possible neighbour Quad in current direction 
        [exist_NQ, refNQ_sec] = refNeighbour(refLeaf_sec,dir);
        if exist_NQ == 1
            % look for neighbour Quad reference in reference array          
            
            idxNQ = findNeighbour(Quadtree,idxLeaf_sec,refLeaf_sec, refNQ_sec);
            if Quadtree.isleaf(idxNQ)
                idxNQa = [idxNQa;idxNQ];
            else
                idxNQchildren = Quadtree.getchildren(idxNQ);
                    if dir == 1 
                      idxNQc1 = idxNQchildren(3);
                      idxNQc2 = idxNQchildren(4);
                
                    elseif dir == 2  
                      idxNQc1 = idxNQchildren(1);
                      idxNQc2 = idxNQchildren(3); 
                    
                    elseif dir == 3  
                      idxNQc1 = idxNQchildren(1);
                      idxNQc2 = idxNQchildren(2); 
                    elseif dir == 4  
                      idxNQc1 = idxNQchildren(2);
                      idxNQc2 = idxNQchildren(4); 
                                
                    end 
                idxNQa = [idxNQa;idxNQc1;idxNQc2];
            end
        end
    end

    

      %get all section numbers of neighbour quad
      sum_secNQ = 0;
      isecNQ = [];
      
      for j = 1: length(idxNQa)                 
          idxsecNQ = find(idxNQa(j,1) == sections(:,3));
          isecNQ = [isecNQ;idxsecNQ];         
          sum_secNQ = sum_secNQ + length(idxsecNQ);
      end
      secNQ{isec}(1,2) = sum_secNQ;
      secNQ{isec}(1,3:2+length(isecNQ)) = isecNQ';
      
   
      
     % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
      %get the neighbour section 
      isecN = 0;
      if ikvo ~= 0 
          sec_NURBS = find(sections(:,4) ~= 0); %get the number of all sections with NURBS curve
          jj = 1;
          while jj <= length(sec_NURBS) && isecN == 0
              idxsecN = sec_NURBS(jj,1);             
              if idxsecN ~= isec
                  if sections(idxsecN,7:end-1) == rot90(sections(isec,7:end-1),2)                 
                      isecN = idxsecN; %section number of the neighbour section 
                  end 
              end
              jj = jj + 1;
          end
          
          
      else %(ikvo == 0)
          
          jj = 1;
          while jj <= length(secNQ{isec}(1,3:end)) && isecN == 0          
              idxsecNQ = secNQ{isec}(1,2+jj);
              if length(sections(idxsecNQ,7:end-1)) == length(sections(isec,7:end-1))                           
                  if sections(idxsecNQ,7:end-2) == rot90(sections(isec,7:end-2),2)            
                      isecN = sections(idxsecNQ,1); %section number of the neighbour section 
                  end 
              end 
              jj = jj + 1;
          end
      end
      secN = [secN;isec,isecN];       
end



%% Elevate the p-/q-Grad

% Get NURBS curve
data = Quadtree.Node{1,1};
NURBS = data{3};
% compute point of the NURBS curve
NURBS_pts = CalculateNURBS(NURBS);

% Extend sections matrix
if ep > 2
    sections = [sections,zeros(length(sections(:,1)),ep-1)];
end

% Insert nodes in P-direction

numSec= size(sections,1);

for isec = 1: numSec             
    ninode = ep - 1; %number of inserted nodes
    ikv = sections(isec,4); %knot vector
    
    if ninode == 0
        continue
    end
    
    coor_StrucNodes = []; %coor_matrix contains the information of all nodes including scalling center for each chosed section
    % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
    if ikv == 0
        %check if this section has been treated in secN
        if sections(isec,6) == 3 %no

            StrucNodes = sections(isec,7:9); %[node_1,node_2,scaling center]
            %extract nodes coordinates of each selected section
            for ii = 1:length(StrucNodes)
                number = StrucNodes(1,ii);
                coor_StrucNodes0 = coor(number,:);
                coor_StrucNodes = [coor_StrucNodes;coor_StrucNodes0];
            end
            %the coordinate of the inserted node in p-direction           
            coor_nodes_p = getLocation(coor_StrucNodes(1,2),coor_StrucNodes(2,2),coor_StrucNodes(1,3),coor_StrucNodes(2,3),ep);
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 0;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            
            % Check if current inserted node is inside the region enclosed by the NURBS curve            
            for j = 1: ninode
                pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor_nodes_p(j,1:2));
                coor(nnode+j,7) = pointInPoly;
            end
                 
            coor(nnode+1:nnode+ninode,2:3) =[coor_nodes_p];
            % update the section matrix 
            nsec_new = 3 + ninode;
            sections(isec,6) = nsec_new;
            sections(isec,6+nsec_new) = sections(isec,9);%shift of the scaling center
            sections(isec,5+nsec_new) = sections(isec,8);%shift of the structral points
            sections(isec,8:7+ninode) = (nnode+1:nnode+ninode);

            %update the ord matrix       
            ord(isec,2) = ep;

            idx_isecN = find(secN(:,1) == isec);
            isecN = secN(idx_isecN,2);
            if isecN ~= 0
                nsecN_new = 3 + ninode;
                sections(isecN,6) = nsecN_new;                 
                sections(isecN,6+nsecN_new) = sections(isecN,9);%shift of the scaling center
                sections(isecN,5+nsecN_new) = sections(isecN,8);%shift of the structural point
                sections(isecN,8:7+ninode) = rot90(sections(isec,8:7+ninode),2);
                ord(isecN,2) = ep;
            end
            nnode = nnode + ninode;
        end        
    end
end
%
%
%
%
%
%
%
%
%
% Extend sections matrix
if eq > 2
    maxpgrad = max(ord(:,2));
    sections = [sections,zeros( size(sections,1), (maxpgrad+1)*(eq-1) )];
end


% Insert nodes for q-direction


numPoly = size(polyElmts,1);
for ipoly = 1: numPoly
    elmt = connectivity{ipoly}(1,6:end); % element connectivity matrix
    kvno = connectivity{ipoly}(2); %knot vector for element
    numSecPoly = polyElmts(ipoly,3);
    secElmts = polyElmts(ipoly,4:3 + numSecPoly);
    ninode = eq - 1; %number of inserted nodes
    
    if ninode == 0
        continue
    end
       
    
    if kvno == 0 
       
        
        %insert nodes in shared edges between sections per poly element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = []; %vertices of this poly element
        numNodesPoly = length(elmt) - 1;
        for inelmt = 1:numNodesPoly
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 0;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            
            % Check if current inserted node is inside the region enclosed by the NURBS curve            
            for j = 1: ninode
                pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor_inodes_q(j,1:2));
                coor(nnode+j,7) = pointInPoly;
            end 
            
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end
        
     
               
        %insert nodes per section within element  
             
        for ise = 1:numSecPoly
            idxse = secElmts(ise);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = []; %nodes of each section except vertices of element
            %insert nodes 
            for indsec = 1:nsec-1                        
                inode = sections(idxse,6+indsec);
                if ismember(inode,nodeElmt(:,1)) == 0 % vertices of element are not included
                    xcoor1 = coor(inode,2);
                    ycoor1 = coor(inode,3);
                    coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                    %update coor matrix
                    coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                    coor(nnode+1:nnode+ninode,4) = 0;  %weight
                    coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                    coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                    
                    % Check if current inserted node is inside the region enclosed by the NURBS curve            
                    for j = 1: ninode
                        pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor_inodes_q(j,1:2));
                        coor(nnode+j,7) = pointInPoly;
                    end 
 
                    coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                    nodeSec0 =[inode,nnode+1:nnode + ninode];
                    nodeSec = [nodeSec;nodeSec0];
                    nnode = nnode + ninode; 
                end
            end
                      
            
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center
            
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
        
    else %kvno ~= 0
        
        %insert nodes in shared edges between sections per element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = [];
        for inelmt = 1:length(elmt)-1
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 0;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            
            % Check if current inserted node is inside the region enclosed by the NURBS curve            
            for j = 1: ninode
                pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor_inodes_q(j,1:2));
                coor(nnode+j,7) = pointInPoly;
            end 
             
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end    
                
  
        %insert nodes per section within element       
        
        for ise = 1:numSecPoly
            idxse = secElmts(ise); %section number within element
            ikv = sections(idxse,4);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = [];
            %insert nodes 
            if ikv == 0                
                for indsec = 1:nsec-1                        
                    inode = sections(idxse,6+indsec);
                    if ismember(inode,nodeElmt(:,1)) == 0
                        xcoor1 = coor(inode,2);
                        ycoor1 = coor(inode,3);
                        coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                        %update coor matrix
                        coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                        coor(nnode+1:nnode+ninode,4) = 0;  %weight
                        coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                        coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                        
                        % Check if current inserted node is inside the region enclosed by the NURBS curve            
                        for j = 1: ninode
                            pointInPoly = isPointInPolygon(NURBS_pts(1:end-1,1:2), coor_inodes_q(j,1:2));
                            coor(nnode+j,7) = pointInPoly;
                        end 
 
                        coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                        nodeSec0 =[inode,nnode+1:nnode + ninode];
                        nodeSec = [nodeSec;nodeSec0];
                        nnode = nnode + ninode; 
                    end
                end 
            end
            
         
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center                
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
    end
end

% adjust the size of sections matrix
maxnsec = max(sections(:,6));
sections(:,6+maxnsec+1:end) = []; 

end