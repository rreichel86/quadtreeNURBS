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
f_plotPolyElmtCurvedEdges = 0; % Plot polygonal elements curve edges
f_splittElmtIntoSec = 1; % Splitt polygonal elements into section

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

%% Plot NURBS curve
figure(1)
hold on
if f_plotNURBS == 1
    plot(NURBS(:,1),NURBS(:,2),'r','LineWidth',3);
    hold on
    plot(controlPoints(1, :), controlPoints(2, :), 'ro','LineWidth',3);
    plot(controlPoints(1, :), controlPoints(2, :), '--','LineWidth',0.5);
end
box on
%% Quadtree decomposition
[Quadtree] = nurbs_brep_quadtree(degree,knots,controlPoints,weights,Boundary);
%% Plots the NURBS contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree,ax)
end
%% Extract polygonal elements 
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots] = extractElements(Quadtree);
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


   
%% Splitt polygonal elements into section

% apply funtion to the 5. and 2. entry of knotVectors cell array. 
% This entries correspond to the degree of the curve and the
% number of knots, respectively.
% ncpoints = (nknot-1) - degree
ncpoints = cellfun(@(x) (x(5)-1)-x(2), knotVectors(1:numKnotVectors));

% compute number of sections and max number of node on any section 
numsec = 0; % number of section 
maxnsec = 0; % max number of nodes on any section
for ielno = 1:numel
    ikv = connectivity{ielno}(2);
    nel = connectivity{ielno}(3);
    if ikv ~= 0
       ncp = ncpoints(ikv); 
       numsec = numsec + nel - ncp + 2;
       nsec = ncp + 1;
    else
       numsec = numsec + nel;
       nsec = 3;
    end
    maxnsec = max(maxnsec,nsec);
end

% sections = [isec, ikv, nsec, node_1,...,node_nsec]
sections = zeros(numsec, maxnsec + 3);
% knots = [ikv, nknots, knot_1,...,knot_nknots]
knots = zeros(numKnotVectors, maxnknots + 2);
% ord = [isec, ikv, pgrad, qgrad]
ord = zeros(numsec, 4);


isec = 0;
if f_splittElmtIntoSec == 1
    for ielno = 1:numel
        kvno = connectivity{ielno}(2);
        nel = connectivity{ielno}(3);
        elmt = connectivity{ielno}(4:end);
        ecoor = coor( elmt(1:end), 2:3);
        wg = coor( elmt(1:end), 4);
        if kvno ~= 0
            nKnot = knotVectors{kvno}(5);
            pgrad = knotVectors{kvno}(2);
            icp = find( elmt == knotVectors{kvno}(3) );
            ecp = find( elmt == knotVectors{kvno}(4) );
            
            knots(kvno,1) = kvno;
            knots(kvno,3) = nKnot;
            knots(kvno,3:2+nKnot) = knotVectors{kvno}(6:end);
            
            ncp = (nKnot - 1) - pgrad;
            
            icp_coor = ecoor(icp,:);
            ecp_coor = ecoor(ecp,:);
            sc_coor = ecoor(end,:);
            % determine orientation
            % ori =  1 for CCW
            % ori = -1 for CW
            ori = orientation(icp_coor,ecp_coor,sc_coor);
            % swap icp and ecp
            if ori == -1
                temp = icp;
                icp = ecp;
                ecp = temp;
            end            
        end
        % determine if a polygonal element is
        % inside or outside the region enclosed 
        % by the NURBS curve
        if sum(coor( elmt, 7)) < 0 % outside
          color = 'red';
        else % inside
          color = 'blue';
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
                
                isec = isec + 1;
                sections(isec,1) = isec;
                sections(isec,2) = 0;
                sections(isec,3) = 3;
                sections(isec,4:6) = elmt(idx);
                
                ord(isec,:) = [isec,0,1,1];
                
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
            end

        % plot polygon that have curve edges
        elseif kvno ~= 0 && icp < ecp && icp == 1
            ii = 1;
            while ii <= nel
                if (ii == ecp) && (icp == 1)
                    a = ii;
                    b = 1;
                    ii = ii + ncp-2;
                    idx = [a:nel,1,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = kvno;
                    sections(isec,3) = ncp + 1;
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,kvno,pgrad,1];
                
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        elseif kvno ~= 0 && icp < ecp && icp ~= 1
            ii = 1;
            while ii <= nel
                if ii == icp
                    a = ii;
                    b = ii + ncp-1;
                    ii = ii + ncp-2;
                    idx = [a:b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = kvno;
                    sections(isec,3) = ncp + 1;
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,kvno,pgrad,1];
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        elseif kvno ~= 0 && icp > ecp && ecp == 1
            ii = 1;
            while ii <= nel
                if (ii == icp) && (ecp == 1)
                    a = ii;
                    b = 1;
                    ii = ii + ncp-2;
                    idx = [a:nel,1,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = ncp + 1;
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,kvno,pgrad,1];
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end

        elseif kvno ~= 0 && icp > ecp && ecp ~= 1
            ii = 1;
            while ii <= nel
                if ii == ecp
                    a = ii;
                    b = ii +  ncp-1;
                    ii = ii + ncp-2;
                    idx = [a:b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = kvno;
                    sections(isec,3) = ncp + 1;
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,kvno,pgrad,1];
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec;
                    sections(isec,2) = 0;
                    sections(isec,3) = 3;
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,0,1,1];
                    
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
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
        kvno = connectivity{ielno}(2);
        nel = connectivity{ielno}(3);
        elmt = connectivity{ielno}(4:3+nel);
        scno = connectivity{ielno}(end);
        ecoor = coor( elmt(1:end), 2:3);
        sc_coor = coor( scno, 2:3);
        wg = coor( elmt(1:end), 4);
        if kvno ~= 0
            nKnot = knotVectors{kvno}(5);
            pgrad = knotVectors{kvno}(2);
            knotVector = knotVectors{kvno}(6:end);
            icp = find( elmt == knotVectors{kvno}(3) );
            ecp = find( elmt == knotVectors{kvno}(4) );
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
        elseif kvno ~= 0 && icp < ecp && icp == 1
            ii = 1;
            while ii <= nel
                if (ii == ecp) && (icp == 1)
                    a = ii;
                    b = 1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor([a:nel,1],:)',wg([a:nel,1])');
                    plot(NURBS_sec(:,1),NURBS_sec(:,2),'b','LineWidth',3);
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
            
        elseif kvno ~= 0 && icp < ecp && icp ~= 1
            ii = 1;
            while ii <= nel
                if ii == icp
                    a = ii;
                    b = ii + (nKnot-1)-pgrad-1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor(a:b,:)',wg(a:b)');
                    plot(NURBS_sec(:,1),NURBS_sec(:,2),'r','LineWidth',3);
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
            
        elseif kvno ~= 0 && icp > ecp && ecp == 1
            ii = 1;
            while ii <= nel
                if (ii == icp) && (ecp == 1)
                    a = ii;
                    b = 1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor([a:nel,1],:)',wg([a:nel,1])');
                    plot(NURBS_sec(:,1),NURBS_sec(:,2),'b','LineWidth',3);
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
            
        elseif kvno ~= 0 && icp > ecp && ecp ~= 1
            ii = 1;
            while ii <= nel
                if ii == ecp
                    a = ii;
                    b = ii +  (nKnot-1)-pgrad-1;
                    ii = ii + (nKnot-1)-pgrad-2;
                    NURBS_sec=CalculateNURBS(pgrad,knotVector,ecoor(a:b,:)',wg(a:b)');
                    plot(NURBS_sec(:,1),NURBS_sec(:,2),'b','LineWidth',3);
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


