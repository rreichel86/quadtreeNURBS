function [Quadtree,nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = refine_quadtree_mesh_2(Quadtree,seedingPoints)
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
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ikv, region, nsec, node_1,...,node_nsec]
% 
% ord ------------------------- section polynomial order
% ord = [isec, pgrad, qgrad]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]

%% Quadtree decomposition
% Get NURBS curve
data = Quadtree.Node{1,1};
NURBS.degree = data{3};
NURBS.knots  = data{4};
NURBS.controlPoints = data{5};
NURBS.weights = data{6};

[Quadtree] = QuadtreeSplit(Quadtree,NURBS,seedingPoints);

[Quadtree] = QuadtreeBalance(Quadtree,NURBS);

%% Extract polygonal elements 
[nnode,coor,numel,connectivity,~,...
 numKnotVectors,knotVectors,idxControlPoints] = extractElements_2(Quadtree);

%% Splitt polygonal elements into section

[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = splittIntoSections_2(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,idxControlPoints);

end
