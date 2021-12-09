function [Quadtree,nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts] = refine_quadtree_mesh(Quadtree,seedingPoints_splitt,seedingPoints_merge,ElevateOrderNumber)
% refine_quadtree_mesh: refine given quadtree based mesh
%
% INPUT: 
% Quadtree ------------------- Quadtree data structure 
% seedingPoints -------------- list of seeding points
%
% OUTPUT:
% Quadtree ------------------- updated Quadtree data structure 
% nnode ---------------------- number of coordinates = number of nodes
% coor ----------------------- nodes coordinates and weights 
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary
%                                             1 - inside 
%                                            -1 - outside
%
% numsec ---------------------- number of sections
% maxnsec --------------------- maximum number of nodes on any section
% sections -------------------- sectionsconnectivity matrix as nsec-tupel of 
%                               nodes, where the first three entries
%                               isec - section number
%                               ipoly - polygonal element number
%                               idxLeaf - index of leaf
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
% 
% ord ------------------------- section polynomial order
% ord = [isec, pgrad, qgrad]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]
%
% polyElmts -------------------- relate sections and polygonal elements
% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
%
% -------------------------------------------------------------------------

figure 
hold on 

%%extract parameters from seedingPoints for different case
if ~exist('ElevateOrderNumber')
    refSeedingPoints_splitt = seedingPoints_splitt;
    ElevateOrderNumber = 0;
else
    [QuadLeaf_splitt,QuadLeaf_merge,refSeedingPoints_splitt,refSeedingPoints_merge] = seedingPoints_settle(seedingPoints_splitt,seedingPoints_merge);
end

%% Quadtree decomposition
% Get NURBS curve
data = Quadtree.Node{1,1};
NURBS = data{3};

[Quadtree] = QuadtreeSplit(Quadtree,NURBS,refSeedingPoints_splitt);

[Quadtree] = QuadtreeBalance(Quadtree,NURBS);
[Quadtree] = check_leaf(Quadtree);

%% Extract polygonal elements 
[nnode,coor,numel,connectivity,~,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);

%% Splitt polygonal elements into section

[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts] = splittIntoSections(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,maxnknots,idxControlPoints);

%% Remain the order from last calculation

if ElevateOrderNumber ~= 0

    [coor,maxnsec,nnode,sections,ord] = remainGrad(Quadtree,nnode,coor,sections,ord,polyElmts,connectivity,QuadLeaf_splitt,QuadLeaf_merge);
                                                  
end

end
