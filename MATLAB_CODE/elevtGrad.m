function [coor,maxnsec,nnode,sections,ord]=elevtGrad(Quadtree,nnode,coor,sections,ord,polyElmts,connectivity,seedingPoints_splitt,seedingPoints_merge)
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


%% get the neighbours of unqualified sections

%seedingPoints_splitt = [isec, isec0, idxLeaf, xcoor, ycoor,pgrad,qgrad]

numSeedingPoints_splitt = length(seedingPoints_splitt(:,1));


secNQ_splitt = cell(numSeedingPoints_splitt,1); %array for all section numbers of neighbour quad
%secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

secN_splitt = []; %array for section numbers of the neighbour section
%secN = [idx_sec,isecN]


for isp = 1: numSeedingPoints_splitt
    idx_sec = seedingPoints_splitt(isp,2); %(old) number of unqualified sections 
    ikvo = sections(idx_sec,4);
    
    idxLeaf_sec = seedingPoints_splitt(isp, 3); %idxLeaf of unqualified sections
      
    refLeaf_sec = Quadtree.Node{idxLeaf_sec,1}{2,1}(1:end);
    secNQ_splitt{isp}(1,1) = isp;
   
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
      secNQ_splitt{isp}(1,2) = sum_secNQ;
      secNQ_splitt{isp}(1,3:2+length(isecNQ)) = isecNQ';
      
   
      
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
          while jj <= length(secNQ_splitt{isp}(1,3:end)) && isecN == 0          
              idxsecNQ = secNQ_splitt{isp}(1,2+jj);
              if length(sections(idxsecNQ,7:end-1)) == length(sections(idx_sec,7:end-1))                           
                  if sections(idxsecNQ,7:end-2) == rot90(sections(idx_sec,7:end-2),2)            
                      isecN = sections(idxsecNQ,1); %section number of the neighbour section 
                  end 
              end 
              jj = jj + 1;
          end
      end
      secN_splitt = [secN_splitt;idx_sec,isecN];       
end

%% get the neighbours of qualified sections
if isempty(seedingPoints_merge) == 0
    
    %seedingPoints_merge = [isec, isec0, idxLeaf, xcoor, ycoor,pgrad,qgrad]
    numSeedingPoints_merge = length(seedingPoints_merge(:,1));

    %array for all section numbers of neighbour quad
    secNQ_merge = cell(numSeedingPoints_merge,1);  %secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

    %array for section numbers of the neighbour section
    secN_merge = [];  %secN = [idx_sec,isecN]


    for isp = 1: numSeedingPoints_merge
        idx_sec = seedingPoints_merge(isp,2); %(old) number of unqualified sections 
        ikvo = sections(idx_sec,4);

        idxLeaf_sec = seedingPoints_merge(isp, 3); %idxLeaf of unqualified sections

        refLeaf_sec = Quadtree.Node{idxLeaf_sec,1}{2,1}(1:end);
        secNQ_merge{isp}(1,1) = isp;

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
          secNQ_merge{isp}(1,2) = sum_secNQ;
          secNQ_merge{isp}(1,3:2+length(isecNQ)) = isecNQ';



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
              while jj <= length(secNQ_merge{isp}(1,3:end)) && isecN == 0          
                  idxsecNQ = secNQ_merge{isp}(1,2+jj);
                  if length(sections(idxsecNQ,7:end-1)) == length(sections(idx_sec,7:end-1))                           
                      if sections(idxsecNQ,7:end-2) == rot90(sections(idx_sec,7:end-2),2)            
                          isecN = sections(idxsecNQ,1); %section number of the neighbour section 
                      end 
                  end 
                  jj = jj + 1;
              end
          end
          secN_merge = [secN_merge;idx_sec,isecN];       
    end
    
else
    secN_merge = [];

end


%% elevate the p-/q-Grad


[coor,nnode,sections,ord] = NodesInsertP(nnode,coor,sections,ord,seedingPoints_splitt,secN_splitt,seedingPoints_merge,secN_merge);

[coor,nnode,sections,ord] = NodesInsertQ(nnode,coor,sections,ord,polyElmts,connectivity,seedingPoints_splitt,seedingPoints_merge);

maxnsec = max(sections(:,6));

end