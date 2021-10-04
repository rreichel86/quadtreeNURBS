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
k_min = 2;
tic 
[Quadtree] = nurbs_brep_quadtree(k_min,NURBS,Boundary);
toc

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
%% Plot sections
figure(3)
hold on
for ie = 1:numsec
    
    kv = sections(ie,2); 
    nsec = sections(ie,4);
    secnd  = sections(ie,5:4+nsec); 
    secx = coor(secnd,2);
    secy = coor(secnd,3);
    
    if kv == 0
        
        patch(secx,secy,'w','FaceAlpha',0,'LineStyle','-','LineWidth',1);
       
    else 
        
        w = knots(kv,2); % weights
        nknots = knots(kv,3); % number of knots
        n = wgt(w,2) - 1; % number of control points or
        %           weights or
        %           NURBS basis functions
        degree = nknots - 1 - (n + 1); % degree of the NURBS curve
        iKnot = knots(kv,4);
        jKnot = knots(kv,5);
        
        controlPoints = zeros(2,n+1);
        controlPoints(1,:) = secx(1:end-1);
        controlPoints(2,:) = secy(1:end-1);
            
        % Compute discrete points of NURBS curve
        NURBS = CalculateNURBS_2(ord(ie,2),iKnot,jKnot,knots(kv,6:end),...
            controlPoints,wgt(w,3:3+n));
            
        % Plot NURBS curve
        plot(NURBS(:,1),NURBS(:,2),'r','LineWidth',2);
        % Plot control points
        plot(controlPoints(1,:), controlPoints(2,:),'b-.');

        secx = [NURBS(:,1); secx(end)];
        secy = [NURBS(:,2); secy(end)];
        patch(secx,secy,'w','FaceAlpha',0,'LineStyle','-');
        
    end    
end                                                                                                                                                 
hold off

%% sections quality 
figure(4)
hold on
colormap(parula(100))
for ie = 1:numsec
    
    kv = sections(ie,2); 
    nsec = sections(ie,4);
    secnd  = sections(ie,5:4+nsec); 
    secx = coor(secnd,2);
    secy = coor(secnd,3);
    
    A = [secx(end); secy(end)];
    C = [secx(end-1); secy(end-1)];
    B = [secx(1); secy(1)];
    
    alpha = qualityTriangle(A,B,C);
    
    patch('Faces',[1:nsec],'Vertices',[secx secy] ,...
        'FaceVertexCData',zeros(nsec,1) + alpha,...
        'FaceColor','interp',...
        'EdgeColor', 'k')
end                                                                                                                                                 
colorbar
hold off