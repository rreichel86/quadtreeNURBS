function [coor,maxnsec,nnode,sections,ord]=elevtGrad(Quadtree,ep,eq,nnode,coor,sections,ord,polyElmts,connectivity,seedingPoints_splitt,seedingPoints_merge)
% elevate order in p-direction and/or q-direction
%
%INPUT:

% Quadtree
% ep = elevated pgrad 
% eq = elevated qgrad 
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
%                               nodes, where the first three entries
%                               isec - section number
%                               ipoly - polygonal element number
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
%
%
% seedingPoints_splitt = [isec,isec0,idxLeaf,xcoor,ycoor]
%                         
%                       isec  --- new section number of the unqualified section
%                       isec0 --- old section number lof the unqualified section 
%                       idxLeaf -- index of Leaf
%                       x/y_coor - x/y coordinate of the scaling center
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


%% get the neighbours of unqualified sections

%seedingPoints = [isec, isec0, idxLeaf, xcoor, ycoor]

num_seedingPoints = length(seedingPoints_splitt(:,1));


secNQ = cell(num_seedingPoints,1); %array for all section numbers of neighbour quad
%secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

secN = []; %array for section numbers of the neighbour section
%secN = [idx_sec,isecN]


for isp = 1: num_seedingPoints
    idx_sec = seedingPoints_splitt(isp,2); %(old) number of unqualified sections 
    ikvo = sections(idx_sec,4);
    
    idxLeaf_sec = seedingPoints_splitt(isp, 3); %idxLeaf of unqualified sections
      
    refLeaf_sec = Quadtree.Node{idxLeaf_sec,1}{2,1}(1:end);
    secNQ{isp}(1,1) = isp;
   
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

    

    % sections = [isec, ipoly, idxLeaf, ikv,region, nsec, node_1,...,node_nsec]    
    % polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]

      %get all section numbers of neighbour quad
      sum_secNQ = 0;
      isecNQ = [];
      
      for j = 1: length(idxNQa)                 
          idxsecNQ = find(idxNQa(j,1) == sections(:,3));
          isecNQ = [isecNQ;idxsecNQ];         
          sum_secNQ = sum_secNQ + length(idxsecNQ);
      end
      secNQ{isp}(1,2) = sum_secNQ;
      secNQ{isp}(1,3:2+length(isecNQ)) = isecNQ';
      
   
      
     % sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
      %get the neighbour section 
      isecN = 0;
      if ikvo ~= 0 
          sec_NURBS = find(sections(:,4) ~= 0); %get the number of all sections with NURBS curve
          jj = 1;
          while jj <= length(sec_NURBS) && isecN == 0
              idxsecN = sec_NURBS(jj,1);             
              if idxsecN ~= idx_sec
                  if sections(idxsecN,7:end-1) == rot90(sections(idx_sec,7:end-1),2)                 
                      isecN = idxsecN; %section number of the neighbour section 
                  end 
              end
              jj = jj + 1;
          end
          
          
      else %(ikvo == 0)
          
          jj = 1;
          while jj <= length(secNQ{isp}(1,3:end)) && isecN == 0          
              idxsecNQ = secNQ{isp}(1,2+jj);
              if length(sections(idxsecNQ,7:end-1)) == length(sections(idx_sec,7:end-1))                           
                  if sections(idxsecNQ,7:end-2) == rot90(sections(idx_sec,7:end-2),2)            
                      isecN = sections(idxsecNQ,1); %section number of the neighbour section 
                  end 
              end 
              jj = jj + 1;
          end
      end
      secN = [secN;idx_sec,isecN];       
end




%% elevate the p-/q-Grad

if ep >= 2 %in p-direction
    [coor,nnode,sections,ord] = NodesInsertP(ep,nnode,coor,sections,seedingPoints_splitt,ord,secN);
end

if eq >= 2 %in q-diretion
    [coor,nnode,sections,ord] = NodesInsertQ(eq,nnode,coor,sections,seedingPoints_splitt,ord,polyElmts,connectivity);
end


if ep >= 2 && eq == 1
    maxnsec = ep + 2;
elseif ep == 1 && eq >= 2
    maxnsec = 3*eq + 1;
elseif ep>=2 && eq>=2
    maxnsec = (ep+1)*eq+1;
end


end