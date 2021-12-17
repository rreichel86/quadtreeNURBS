function [coor,maxnsec,nnode,sections,ord]=remainGrad(Quadtree,nnode,coor,sections,ord,polyElmts,connectivity,QuadLeaf_splitt,QuadLeaf_merge)
%remain the p-/qgrad of quad from last calculation
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
%                              nodes, where the first five entries
%                              iel - element number
%                              ikv - knot vector number
%                              idxLeaf - index of Leaf
%                              which_region - region number
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]

% QuadLeaf = [iQuad,quadkv,maxpgrad,maxqgrad]
%             iQuad: Leaf number of (un)qualified quad;
%             quadkv: 0(no knot vector),1(knot vector);
%             maxpgrad: max. pgrad from all sections of current quad
%             maxqgrad: max. qgrad from all sections of current quad
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




%% new seedingPoints based on original type

%number of unqualified quad
nQuadLeaf_splitt = size(QuadLeaf_splitt,1);
newSeedingPoints_splitt = [];
for i = 1: nQuadLeaf_splitt
    iQuadLeaf = QuadLeaf_splitt(i,1);
    pgrad = QuadLeaf_splitt(i,3);
    qgrad = QuadLeaf_splitt(i,4);
    
    %check if current quad has children
    idxChildren = Quadtree.getchildren(iQuadLeaf);
    numChildren = length(idxChildren);
    if numChildren ~= 0 %yes
        for ii = 1: numChildren
            ichild = idxChildren(ii);
            idxsecs = find(sections(:,3)==ichild);
            for iii = 1:length(idxsecs)
                isec = idxsecs(iii);
                nsec = sections(isec,6);
                ikv = sections(isec,4);
                sc = sections(isec,6+nsec);
                sc_coor = coor(sc,2:3);
                newSeedingPoints_splitt = [newSeedingPoints_splitt;isec,isec,ichild,ikv,sc_coor,pgrad-1,qgrad-1];
            end
        end
        
    else %numChildren == 0, no     %p/q-grad of all sections in current quad will be unified into one ???
        idxsecs = find(sections(:,3)==iQuadLeaf);        
        for ii = 1:length(idxsecs)
            isec = idxsecs(ii);
            nsec = sections(isec,6);
            ikv = sections(isec,4);
            sc = sections(isec,6+nsec);
            sc_coor = coor(sc,2:3);
            newSeedingPoints_splitt = [newSeedingPoints_splitt;isec,isec,iQuadLeaf,ikv,sc_coor,pgrad-1,qgrad-1];
        end       
    end
end

newSeedingPoints_merge = [];

nQuadLeaf_merge = size(QuadLeaf_merge,1);
for i = 1: nQuadLeaf_merge
    iQuadLeaf = QuadLeaf_merge(i,1);
    pgrad = QuadLeaf_merge(i,3);
    qgrad = QuadLeaf_merge(i,4);
    
   %check if current quad has already be treated in QuadLeaf_splitt
   if ismember(iQuadLeaf,QuadLeaf_splitt(:,1)) == 1 %yes
       continue
   end

    %check if current quad has children
    idxChildren = Quadtree.getchildren(iQuadLeaf);
    numChildren = length(idxChildren);
    if numChildren ~= 0 %yes
        for ii = 1: numChildren
            ichild = idxChildren(ii);
            idxsecs = find(sections(:,3)==ichild);
            for iii = 1:length(idxsecs)
                isec = idxsecs(iii);
                nsec = sections(isec,6);
                ikv = sections(isec,4);
                sc = sections(isec,6+nsec);
                sc_coor = coor(sc,2:3);
                newSeedingPoints_merge = [newSeedingPoints_merge;isec,isec,ichild,ikv,sc_coor,pgrad,qgrad];
            end
        end

    else %numChildren == 0, no    %p/q-grad of all sections in current quad will be unified into one ???      
        idxsecs = find(sections(:,3)==iQuadLeaf);        
        for ii = 1:length(idxsecs)
            isec = idxsecs(ii);
            nsec = sections(isec,6);
            ikv = sections(isec,4);
            sc = sections(isec,6+nsec);
            sc_coor = coor(sc,2:3);
            newSeedingPoints_merge = [newSeedingPoints_merge;isec,isec,iQuadLeaf,ikv,sc_coor,pgrad,qgrad];
        end       
    end
end


%% get the neighbours of unqualified sections

numSeedingPoints_splitt = size(newSeedingPoints_splitt,1);

secNQ_splitt = cell(numSeedingPoints_splitt,1); %array for all section numbers of neighbour quad
%secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

secN_splitt = []; %array for section numbers of the neighbour section
%secN = [idx_sec,isecN]


for isp = 1: numSeedingPoints_splitt
    idx_sec = newSeedingPoints_splitt(isp,2); %(old) number of unqualified sections 
    ikvo = sections(idx_sec,4);
    
    idxLeaf_sec = newSeedingPoints_splitt(isp, 3); %idxLeaf of unqualified sections
      
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

    
%seedingPoints_merge = [isec, isec0, idxLeaf, xcoor, ycoor,pgrad,qgrad]
numSeedingPoints_merge = size(newSeedingPoints_merge,1);

%array for all section numbers of neighbour quad
secNQ_merge = cell(numSeedingPoints_merge,1);  %secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

%array for section numbers of the neighbour section
secN_merge = [];  %secN = [idx_sec,isecN]


for isp = 1: numSeedingPoints_merge
    idx_sec = newSeedingPoints_merge(isp,2); %(old) number of unqualified sections 
    ikvo = sections(idx_sec,4);

    idxLeaf_sec = newSeedingPoints_merge(isp, 3); %idxLeaf of unqualified sections

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
    




%% elevate the p-/q-Grad


[coor,nnode,sections,ord] = NodesInsertP(nnode,coor,sections,ord,newSeedingPoints_splitt,secN_splitt,newSeedingPoints_merge,secN_merge);

[coor,nnode,sections,ord] = NodesInsertQ(nnode,coor,sections,ord,polyElmts,connectivity,newSeedingPoints_splitt,newSeedingPoints_merge);

maxnsec = max(sections(:,6));



end