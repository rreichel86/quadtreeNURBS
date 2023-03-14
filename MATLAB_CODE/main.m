clearvars;
close all;
clc;
config;

% Select example
example_nro = 1;
%  |      Nr.     |      Description        |
%  |----------------------------------------|
%  |       1      |       Moby-Dick         |
%  |       2      |     Circumference       |
%  |       3      |    Double circumf.      |
%  |       4      |      Flat shape         |
%  |       5      |       Heart             |
%  |       6      |  A quarter circumf.     |

% Plot options
f_plotNURBS = 1; % Plot NURBS curve
f_plotLeaves = 1; % Plot the NURBS contained in each leaf separately

% Initialization
% ==============
% Control points input should be a matrix with dimensions
% (coordinates, nPoints), and should be as follows: first row corresponds to
% the x coordinate, second to the y coordinate. Each control point
% consequently lies in its own column.

% Obtains the selected NURBS definition
[NURBS,Boundary,FileName] = NURBS_parameters(example_nro);

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
    set(gca,'TickLabelInterpreter','latex','FontSize',20,'FontName','Times');
    axis square
    axis([xmin xmax ymin ymax])
    hold on

    % Plot control polygon
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'b-.','LineWidth',1);
    % Plot control points
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'o','Color','red','MarkerFaceColor','r','MarkerSize',6);
    % Plot curve
    plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2.0);
    % Plot boundary
    patch(Boundary(1,:),Boundary(2,:), 'w','FaceAlpha',0,'LineWidth',1.5)

end

% currentFileName = strcat('../Examples/',FileName,'0');
% print(currentFileName,'-dpng');

%% Quadtree decomposition
k_min = 0;
[Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary);

%% Plots the NURBS segments contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree)

%     currentFileName = strcat('../Examples/',FileName,'1');
%     print(currentFileName,'-dpng');
end

%% Extract polygonal elements
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);

%% Splitt polygonal elements into section
[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts] = splittIntoSections(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,maxnknots,idxControlPoints);

