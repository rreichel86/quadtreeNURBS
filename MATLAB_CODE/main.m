clearvars;
close all;
clc;
config;

% Select example
example_nro = 7;
%      Nr.     |      Description        |    Section     
%---------------------------------------------------------
%       1      |       Moby-Dick         |      3         
%       2      |     Circumference       |     4.
%       3      |    Double circumf.      |     4.1
%       4      |      Flat shape         |     4.1
%       5      |  Limit double circum.   |     4.2
%       6      |    Degenerated case     |     4.2
%       7      |       Heart             |      -


% Plot options
f_plotNURBS = 1; % Plot NURBS curve
f_plotLeaves = 0; % Plot the NURBS contained in each leaf separately
f_plotPolyElmt = 0; % Plot polygonal elements 
f_splittElmtIntoSec = 1; % Splitt polygonal elements into section

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

%% Plot NURBS curve
figure(1)
hold on
if f_plotNURBS == 1
    plot(NURBS_pts(:,1),NURBS_pts(:,2),'r','LineWidth',2.5);
    hold on
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :), 'b-.','LineWidth',1);
    plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :), 'o','Color','red','MarkerFaceColor','r','MarkerSize',8);
    hold on 
    patch(Boundary(1,:),Boundary(2,:), 'w','FaceAlpha',0)
end
box on


%% Quadtree decomposition
k_min = 2;
[Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary);

%% Plots the NURBS segments contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree)
end
%% Extract polygonal elements 
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);

%% Plot polygonal elements 
% if f_plotPolyElmt == 1
%     for ielno = 1:numel
%         
%         region_nro = connectivity{ielno}(3);
%         nel = connectivity{ielno}(6);
%         elmt = connectivity{ielno}(7:6+nel);
%         sc = connectivity{ielno}(end);
%         
%         ex = coor( elmt, 2);
%         ey = coor( elmt, 3);
%         
%         sc_x = coor( sc, 2);
%         sc_y = coor( sc, 3);
%         
%         % determine if a polygonal element is
%         % inside or outside the region enclosed 
%         % by the NURBS curve
%         if region_nro == 1
%         %if sum(coor( elmt, 7)) < 0 % outside
%           color = 'red';
%         else % inside
% %           continue
%           color = 'blue';
%         end
%         
%         patch(ex,ey, color,'FaceAlpha',.5)
%         hold on
%         plot(sc_x, sc_y, 'k*')
%     end
% end

%% Splitt polygonal elements into section

[nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts] = splittIntoSections(nnode,coor,numel,connectivity,...
                                                                    numKnotVectors,knotVectors,maxnknots,idxControlPoints);