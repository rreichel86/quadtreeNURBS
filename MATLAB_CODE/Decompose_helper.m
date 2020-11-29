function [Quadtree] = Decompose_helper(Quadtree,NURBS,idxQ)
% Decompose_helper: prepares input for decompose SR. 
% It creates the auxiliar variables needed for calling the decompose SR.

Quad = Quadtree.Node{idxQ,1}{10};
LocQ = ref2loc(Quadtree.Node{idxQ,1}{2});
idxFather = Quadtree.Parent(idxQ);
Q_aux = Quadtree.Node{idxFather,1}{2};
k = length(LocQ);
l = LocQ(1);
pos_aux = LocQ;

% decompose SR
[Quadtree] = decompose(Quadtree,Quad,NURBS,l,k,pos_aux,Q_aux);

end