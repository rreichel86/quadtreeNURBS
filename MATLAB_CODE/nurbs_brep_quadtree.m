function [Quadtree] = nurbs_brep_quadtree(degree,knots,controlPoints,weights,Boundary)


close all;
clc;
config;

% CalculateNURBS computes the curve
NURBS = CalculateNURBS(degree,knots,controlPoints,weights);

% First plotting
figure(1);
hold on;
axis(ax);
plot(NURBS(1,:),NURBS(2,:),'r','LineWidth',3);
plot(controlPoints(1, :), controlPoints(2, :), 'ro','LineWidth',3);
plot(controlPoints(1, :), controlPoints(2, :), '--','LineWidth',0.5);
box on

%Count of control points in the root
nPoints = checkQuad( Boundary(:,2)',  Boundary(:,4)', Boundary(:,2)',...
    controlPoints);

% Split the root if there is more than one point in the root
if nPoints > 1
    %Setting tree
    Quadtree = tree('root');
    data=Quadtree.Node{1,1};
    data={data;[]};
    Quadtree = Quadtree.set(1, data);
    
    l=1; % l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initialize at 1
    k=0; % k: level of decomposition, equal to 0 at the root
    pos_aux=[]; Q_aux=[0]; %auxiliar arrays
    
    [Quadtree] = decompose(Quadtree,Boundary,NURBS ,controlPoints, knots,...
        weights, degree, l, k,pos_aux, Q_aux,Boundary);
end

[Quadtree] = Star_Shape(Quadtree,controlPoints,NURBS, knots, weights,...
    degree, Boundary);


[Quadtree] = Balance_Quadtree(Quadtree,NURBS,controlPoints, knots, ...
    weights, degree,Boundary);


end
