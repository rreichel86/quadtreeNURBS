clear all;
close all;
clc;

% Controls:
%   Figure Nr. |      Description        |    Section     
%---------------------------------------------------------
%       1      |       Moby-Dick         |      3         
%       2      |     Circumference       |     4.1
%       3      |    Double circumf.      |     4.1
%       4      |      Flat shape         |     4.1
%       5      |  Limit double circum.   |     4.2
%       6      |    Degenerated case     |     4.2
fig=3;



% Initialization
% ==============
% Control points input should be a matrix with dimensions
% (coordinates,nPoints), and should be as follows: first row corresponds to
% the x coordinate, second to the y coordinate, . Each control point
% consequently lies in its own column.

% Obtains the selected NURBS definition
[degree,knots,controlPoints,weights,ax,Boundary]=NURBS_parameters(fig);

% CalculateNURBS computes the curve
NURBS=CalculateNURBS(degree,knots,controlPoints,weights);

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

%hold off
                        
[Quadtree]=Star_Shape(Quadtree,controlPoints,NURBS, knots, weights,...
                        degree, Boundary);


[Quadtree]=Balance_Quadtree(Quadtree,NURBS,controlPoints, knots, ...
                            weights, degree,Boundary);

% close all the figures and plots the NURBS contained in each leaf separately
%plot_leaf(Quadtree,ax)
[coor,connectivity,nnode,maxnel,numel,kv_element,kv_num,maxnk]=extract_leaf(Quadtree);

%% pick first element of element_nodes 
for ielno = 1:numel
    elmt = connectivity{ielno,:};
    ex = coor( elmt(3:end), 2);
    ey = coor( elmt(3:end), 3);
    
    patch(ex,ey, 'red','FaceAlpha',.5)
end 
