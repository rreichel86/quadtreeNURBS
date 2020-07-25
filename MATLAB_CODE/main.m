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


% Plot options
f_plotLeaves = 0; % Plots the NURBS contained in each leaf separately
f_plotPolyElmt = 1; % Plot polygonal elements
f_plotPolyElmtCurvedEdges = 0; % Plot polygonal elements curve edges
f_splittElmtIntoSec = 0; % Splitt polygonal elements into section

% Initialization
% ==============
% Control points input should be a matrix with dimensions
% (coordinates,nPoints), and should be as follows: first row corresponds to
% the x coordinate, second to the y coordinate. Each control point
% consequently lies in its own column.

% Obtains the selected NURBS definition
[degree,knots,controlPoints,weights,ax,Boundary]=NURBS_parameters(example_nro);

% compute point of the NURBS curve
NURBS = CalculateNURBS(degree,knots,controlPoints,weights);

% Plot NURBS curve
figure(1)
plot(NURBS(:,1),NURBS(:,2),'r','LineWidth',3);
hold on
plot(controlPoints(1, :), controlPoints(2, :), 'ro','LineWidth',3);
plot(controlPoints(1, :), controlPoints(2, :), '--','LineWidth',0.5);
box on

%% Quadtree decomposition
[Quadtree] = nurbs_brep_quadtree(degree,knots,controlPoints,weights,Boundary);

%% Extract polygonal elements 
[nnode,coor,numel,connectivity,maxnel,kv_element,kv_num,maxnk]=extractElements(Quadtree);
%% Plot Nodes

figure(2)
plot(NURBS(:,1),NURBS(:,2),'b','LineWidth',1);
hold on
plot(coor(1:nnode-numel-1,2),coor(1:nnode-numel-1,3),'.r') ;
hold on;
plot(coor(nnode-numel:nnode,2),coor(nnode-numel:nnode,3),'.k') ;
hold on;
% text(coor(:,2), coor(:,3), num2str(coor(:,1)));

% Compute bounding box that enclosed the NURBS curve
x_min = min(NURBS(:,1));
x_max = max(NURBS(:,1));
y_min = min(NURBS(:,2));
y_max = max(NURBS(:,2));

% Plot bounding box
% plot( [x_min,x_max,x_max,x_min,x_min], [y_min,y_min,y_max,y_max,y_min], '-k')
% hold on;

% loop over nodes, excluding control points
for ii = find(coor(:,5) == 1)'
    
    % Check if current node is inside the bounding box
    if ( isPointInQuad([x_min,y_min], [x_max,y_max], coor(ii,2:3)) == 1 )
        
        
%         plot(coor(ii,2),coor(ii,3),'xc')
%         hold on;
        
        % Check if current node is also inside 
        % the region enclosed by the NURBS curve
        pointInPoly = isPointInPolygon(NURBS(1:end-1,1:2), coor(ii,2:3));
        coor(ii,7) = pointInPoly;
%         if (  pointInPoly > 0 )
%             
%             
%             plot(coor(ii,2),coor(ii,3),'*b')
%             hold on;
%             
%         end
    end
end

%% Plots the NURBS contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree,ax)
end
   
%% Splitt polygonal elements into section
if f_splittElmtIntoSec == 1
    for ielno = 1:numel
        kvno = kv_element{ielno}(2);
        nel = kv_element{ielno}(3);
        elmt = kv_element{ielno}(4:end);
        ecoor = coor( elmt(1:end), 2:3);
        wg = coor( elmt(1:end), 4);
        if kvno ~= 0
            nKnot = kv_num{kvno}(5);
            pgrad = kv_num{kvno}(2);
            knotVector = kv_num{kvno}(6:end);
            iknot = find( elmt == kv_num{kvno}(3) );
            eknot = find( elmt == kv_num{kvno}(4) );
            
            iknot_coor = ecoor(iknot,:);
            eknot_coor = ecoor(eknot,:);
            sc_coor = ecoor(end,:);
            % determine orientation
            % ori =  1 for CCW
            % ori = -1 for CW
            ori = orientation(iknot_coor,eknot_coor,sc_coor);
            
            % swap iknot and eknot
            if ori == -1
                temp = iknot;
                iknot = eknot;
                eknot = temp;
            end
            
            
        end
        
        % plot polygon that dont have curve edges
        if kvno == 0
            for ii = 1:nel
                if ii ~= nel
                    a = ii;
                    b = ii+1;
                else
                    a = nel;
                    b = 1;
                end
                idx = [a,b,nel+1];
                patch(ecoor(idx,1).', ecoor(idx,2).', 'green','FaceAlpha',.5)

                hold on
            end

        % plot polygon that have curve edges
        elseif kvno ~= 0 && iknot < eknot && iknot == 1
            ii = 1;
            while ii <= nel
                if (ii == eknot) && (iknot == 1)
                    a = ii;
                    b = 1;
                    ii = ii + pgrad - 1;
                    idx = [a:nel,1,nel+1];
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', 'blue','FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        elseif kvno ~= 0 && iknot < eknot && iknot ~= 1
            ii = 1;
            while ii <= nel
                if ii == iknot
                    a = ii;
                    b = ii + pgrad;
                    ii = ii + pgrad - 1;
                    idx = [a:b,nel+1];
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', 'blue','FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        elseif kvno ~= 0 && iknot > eknot && eknot == 1
            ii = 1;
            while ii <= nel
                if (ii == iknot) && (eknot == 1)
                    a = ii;
                    b = 1;
                    ii = ii + pgrad - 1;
                    idx = [a:nel,1,nel+1];
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', 'blue','FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end

        elseif kvno ~= 0 && iknot > eknot && eknot ~= 1
            ii = 1;
            while ii <= nel
                if ii == eknot
                    a = ii;
                    b = ii + pgrad;
                    ii = ii + pgrad - 1;
                    idx = [a:b,nel+1];
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', 'blue','FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        end
    end
end

%% Plot polygonal elements 
if f_plotPolyElmt == 1
    for ielno = 1:numel
        nel = connectivity{ielno}(3);
        elmt = connectivity{ielno}(4:3+nel);
        ex = coor( elmt, 2);
        ey = coor( elmt, 3);
        
        % determine if a polygonal element is
        % inside or outside the region enclosed 
        % by the NURBS curve
        if sum(coor( elmt, 7)) < 0 % outside
          color = 'red';
        else % inside
          color = 'blue';
        end
        
        patch(ex,ey, color,'FaceAlpha',.5)
    end
end
%% Plot polygonal elements with curve edges
if f_plotPolyElmtCurvedEdges == 1
    for ielno = 1:numel
        kvno = kv_element{ielno}(2);
        nel = kv_element{ielno}(3);
        elmt = kv_element{ielno}(4:4+nel-1);
        scno = kv_element{ielno}(end);
        ecoor = coor( elmt(1:end), 2:3);
        sc_coor = coor( scno, 2:3);
        wg = coor( elmt(1:end), 4);
        if kvno ~= 0
            nKnot = kv_num{kvno}(5);
            pgrad = kv_num{kvno}(2);
            knotVector = kv_num{kvno}(6:end);
            iknot = find( elmt == kv_num{kvno}(3) );
            eknot = find( elmt == kv_num{kvno}(4) );
        end
        
        % plot polygon that dont have curve edges
        if kvno == 0
            %         for ii = 1:nel
            %             if ii ~= nel
            %                 a = ii;
            %                 b = ii+1;
            %             else
            %                 a = nel;
            %                 b = 1;
            %             end
            %             plot(ecoor([a,b],1).', ecoor([a,b],2).', 'c-')
            %
            %             hold on
            %         end
            
            patch(ecoor(:,1),ecoor(:,2), 'green','FaceAlpha',.5)
            
        % plot polygon that dont have curve edges
        elseif kvno ~= 0 && iknot < eknot && iknot == 1
            ii = 1;
            while ii <= nel
                if (ii == eknot) && (iknot == 1)
                    a = ii;
                    b = 1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor([a:nel,1],:)',wg([a:nel,1])');
                    plot(NURBS_sec(1,:),NURBS_sec(2,:),'b','LineWidth',3);
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                elseif ii == nel
                    a = nel;
                    b = 1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                end
                hold on
                ii = ii + 1;
            end
            
        elseif kvno ~= 0 && iknot < eknot && iknot ~= 1
            ii = 1;
            while ii <= nel
                if ii == iknot
                    a = ii;
                    b = ii + (nKnot-1)-pgrad-1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor(a:b,:)',wg(a:b)');
                    plot(NURBS_sec(1,:),NURBS_sec(2,:),'r','LineWidth',3);
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'r-')
                elseif ii == nel
                    a = nel;
                    b = 1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'r-')
                end
                hold on
                ii = ii + 1;
            end
            
        elseif kvno ~= 0 && iknot > eknot && eknot == 1
            ii = 1;
            while ii <= nel
                if (ii == iknot) && (eknot == 1)
                    a = ii;
                    b = 1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor([a:nel,1],:)',wg([a:nel,1])');
                    plot(NURBS_sec(1,:),NURBS_sec(2,:),'b','LineWidth',3);
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                elseif ii == nel
                    a = nel;
                    b = 1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                end
                hold on
                ii = ii + 1;
            end
            
        elseif kvno ~= 0 && iknot > eknot && eknot ~= 1
            ii = 1;
            while ii <= nel
                if ii == eknot
                    a = ii;
                    b = ii +  (nKnot-1)-pgrad-1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor(a:b,:)',wg(a:b)');
                    plot(NURBS_sec(1,:),NURBS_sec(2,:),'b','LineWidth',3);
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                elseif ii == nel
                    a = nel;
                    b = 1;
                    plot(ecoor([a,b],1).', ecoor([a,b],2).', 'b-')
                end
                hold on
                ii = ii + 1;
            end
            
        end
        
    end
end


