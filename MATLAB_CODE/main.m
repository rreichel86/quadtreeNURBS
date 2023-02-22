clearvars;
close all;
clc;
config;

% Select example
example_nro = 2;
%      Nr.     |      Description        |    Section
%---------------------------------------------------------
%       1      |       Moby-Dick         |      3
%       2      |     Circumference       |     4.
%       3      |    Double circumf.      |     4.1
%       4      |      Flat shape         |     4.1
%       5      |  Limit double circum.   |     4.2
%       6      |    Degenerated case     |     4.2
%       7      |       Heart             |      -
%       8      |  A quarter circumf.     |

% Plot options
f_plotNURBS = 1; % Plot NURBS curve
f_plotLeaves = 1; % Plot the NURBS contained in each leaf separately
f_plotPolyElmt = 0; % Plot polygonal elements
f_splittElmtIntoSec = 0; % Splitt polygonal elements into section

% Initialization
% ==============
% Control points input should be a matrix with dimensions
% (coordinates,nPoints), and should be as follows: first row corresponds to
% the x coordinate, second to the y coordinate. Each control point
% consequently lies in its own column.

% Obtains the selected NURBS definition
[NURBS,Boundary] = NURBS_parameters(example_nro);

% compute point of the NURBS curve
NURBS_pts = CalculateNURBS(NURBS);

numRef = 0;
NURBS = hRefinement1d(NURBS,numRef);

%% Plot NURBS curve
if f_plotNURBS == 1
    
    xmin = min(Boundary(1,:));
    xmax = max(Boundary(1,:));
    ymin = min(Boundary(2,:));
    ymax = max(Boundary(2,:));
    
    figure(1)
    xticks([])
    yticks([])
    daspect([1 1 1])
    box on
    set(gca,'TickLabelInterpreter','latex','FontSize',18,'FontName','Times');
    axis square
    axis([xmin xmax ymin ymax])
    hold on
    
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'b-.','LineWidth',1);
    hold on 
    plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2.0);
    hold on
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'o','Color','red','MarkerFaceColor','r','MarkerSize',6);
    hold on
    patch(Boundary(1,:),Boundary(2,:), 'w','FaceAlpha',0,'LineWidth',1.5)
    
end

%% Quadtree decomposition
k_min = 0;
tic
[Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary);
toc

[Quadtree] = check_leaf(Quadtree);
%% Plots the NURBS segments contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree)
end
%% Extract polygonal elements
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);

%% Splitt polygonal elements into section
[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts] = splittIntoSections(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,maxnknots,idxControlPoints);

