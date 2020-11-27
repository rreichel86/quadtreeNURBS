function [Quadtree,nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = nurbs_quadtree_mesh_2(k_min,NURBS,Boundary)
% nurbs_quadtree_mesh_2: generate a quadtree based mesh 
%
% INPUT: 
% k_min ---------------------- minimum level of decomposition
% NURBS definition
% NURBS.degree --------------- NURBS degree
% NURBS.knots ---------------- NURBS knot vector
% NURBS.controlPoints -------- NURBS control points
% NURBS.weights -------------- NURBS weights
% Boundary ------------------- Outer boundary 
% 
% OUTPUT:
% Quadtree ------------------- Quadtree data structure
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
% sections -------------------- sections connectivity matrix as nsec-tupel of 
%                               nodes, where the first three entries
%                               isec - section number
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ikv, region, nsec, node_1,...,node_nsec]
% 
% ord ------------------------- section polynomial order
% ord = [isec, pgrad, qgrad, pos, Npos]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]
%
% -------------------------------------------------------------------------

%% Quadtree decomposition
[Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary);

%% Extract polygonal elements 
[nnode,coor,numel,connectivity,~,...
 numKnotVectors,knotVectors,idxControlPoints] = extractElements_2(Quadtree);

%% Splitt polygonal elements into section

[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = splittIntoSections_2(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,idxControlPoints);

end
