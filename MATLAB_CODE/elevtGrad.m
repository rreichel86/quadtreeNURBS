function [coor,maxnsec,nnode,sections,ord]=elevtGrad(ep,eq,nnode,coor,sections,seedingPoints_splitt,ord,secN,polyElmts,connectivity)
% elevate the grad in p-direction and/or q-direction
%
%INPUT:
%ep = elevated pgrad 
%eq = elevated qgrad 
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
%                               ikv - knot vector number
%                               iel - element number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, idxLeaf, ikv, iel, region, nsec, node_1,...,node_nsec]
%
%
% seedingPoints_splitt = [isec,isec0,idxLeaf,xcoor,ycoor,c]
%                         
%                       isec  --- new section number of the unqualified section
%                       isec0 --- old section number lof the unqualified section 
%                       idxLeaf -- index of Leaf
%                       x/y_coor - x/y coordinate of the scaling center
%                       c  ----   error_measure
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
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 


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