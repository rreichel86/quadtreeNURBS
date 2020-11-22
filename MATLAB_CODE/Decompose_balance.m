function [Quadtree]=Decompose_balance(Quadtree,NURBS,i,Loc_Current,idx_Father,Boundary)
% Decompose_balance functions prepares input for decompose function. It
% creates the auxiliar variables needed for calling the decompose function

Quad=Quadtree.Node{i,1}{10};
Q_aux=Quadtree.Node{idx_Father,1}{2};
k=length(Loc_Current);
l=Loc_Current(1);
ii=Loc_Current(end);
pos_aux=Loc_Current;



[Quadtree] = decompose(Quadtree,Quad,NURBS,l,k,pos_aux,Q_aux,Boundary);


end
