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

%% Plots the NURBS contained in each leaf separately
if f_plotLeaves == 1
    plot_leaf(Quadtree,ax)
end
%% Extract polygonal elements 
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);
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

% compute number of sections, 
% max number of node on any section and
% max number of control points on any section 
%     with boundary defined by a NURBS
numsec = 0; % number of sections
numsec_w_NURBS = 0; % number of sections 
%                     with boundary defined by a NURBS
maxnsec = 0; % max number of nodes on any section
maxncp = 0; % max number of control points
for ielno = 1:numel
    ikv = connectivity{ielno}(2);
    nel = connectivity{ielno}(3);
    if ikv ~= 0
       ncp = ncpoints(ikv); 
       numsec = numsec + nel - ncp + 2;
       numsec_w_NURBS = numsec_w_NURBS + 1;
       nsec = ncp + 1;
       maxncp = max(maxncp, ncp);
    else
       numsec = numsec + nel;
       nsec = 3;
    end
    maxnsec = max(maxnsec,nsec);
end

% sections = [isec, ikv, nsec, node_1,...,node_nsec]
sections = zeros(numsec, maxnsec + 3);
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
knots = zeros(numsec_w_NURBS, maxnknots + 5);
% ord = [isec, pgrad, qgrad]
ord = zeros(numsec, 3);
% weights = [iw, nweights, weight_1,...,weigth_nweigths]
wgt = zeros(numsec_w_NURBS, maxncp + 2);

isec = 0;
isec_w_NURBS = 0;
if f_splittElmtIntoSec == 1
    
    
    
    
    for ielno = 1:numel
        kvno = connectivity{ielno}(2); % knot vector number
        nel = connectivity{ielno}(3); % number of nodes per element
        elmt = connectivity{ielno}(4:end); % element connectivity matrix
        ecoor = coor( elmt(1:end), 2:3);
        wg = coor( elmt(1:end), 4);
       
        % determine if a polygonal element is
        % inside or outside the region enclosed 
        % by the NURBS curve
        if sum(coor( elmt, 7)) < 0 % outside
          color = 'red';
        else % inside
          color = 'blue';
        end
        
        if kvno ~= 0
            
            
            nKnot = knotVectors{kvno}(5);
            pgrad = knotVectors{kvno}(2);
            icp = find( elmt == knotVectors{kvno}(3) );
            ecp = find( elmt == knotVectors{kvno}(4) );
            
            
            
            % number of control points
            ncp = (nKnot - 1) - pgrad;
            
            % scaling center coordinates
            sc_coor = ecoor(end,:);
            
            % control points indices
            idxCtrlP = idxControlPoints{kvno}(3:end);
            
            poly = zeros(ncp+1,2);
            
            
            poly(1:ncp,:) = coor(idxCtrlP,2:3);
            poly(ncp+1,:) = sc_coor; 
            
            % determine signed polygon area
            % area =  1 for CCW
            % area = -1 for CW
            area = poly(ncp+1,1) * poly(1,2) - poly(1,1) * poly(ncp+1,2);
            
            for i = 1: ncp
                area = area + poly(i,1) * poly(i+1,2);
                area = area - poly(i+1,1) * poly(i,2);
            end
            
            % - area is positive, OP = 1. 
            % vertices are arranged CCW
            % - area is negative, OP = -1. 
            % vertices are arranged CW
            if area > 0
                OP = 1
            else
                OP = -1
            end
            
            % if vertices are arranged CW, 
            % swap icp and ecp 
            if OP == -1
                temp = icp;
                icp = ecp;
                ecp = temp;
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
                
                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = 0; % knot vector number (ikv)
                sections(isec,3) = 3; % number of nodes per section (nsec)
                sections(isec,4:6) = elmt(idx); 
                
                ord(isec,:) = [isec,1,1];
                
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
            end

            
            
        % plot polygon that have curve edges
        elseif kvno ~= 0 && icp < ecp && icp == 1
            ii = 1;
            while ii <= nel
                if ii == icp
                    a = 1;
                    b = ii + ncp-1;
                    ii = ii + ncp-2;
                    idx = [a:b,nel+1];
                    
                    isec = isec + 1;
                    isec_w_NURBS = isec_w_NURBS + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                    sections(isec,3) = ncp + 1; % number of nodes per section (nsec)
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,pgrad,1];
                    
                    % knot vector number (ikv)
                    knots(isec_w_NURBS,1) = isec_w_NURBS;
                    % weights number (iw)
                    knots(isec_w_NURBS,2) = isec_w_NURBS;
                    % number of knots (nknots)
                    knots(isec_w_NURBS,3) = nKnot;
                    % initial knot value (iknot)
                    knots(isec_w_NURBS,4) = knotVectors{kvno}(6);
                    % final knot value (jknot)
                    knots(isec_w_NURBS,5) = knotVectors{kvno}(end);
                    % knots
                    knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                    
                    % weights number (iw)
                    wgt(isec_w_NURBS,1) = isec_w_NURBS;
                    % number of control points / weights (nweights)
                    wgt(isec_w_NURBS,2) = ncp;
                    % weights
                    wgt(isec_w_NURBS,3:2+ncp) = coor(elmt(idx(1:end-1)),4);
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
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
                    isec_w_NURBS = isec_w_NURBS + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                    sections(isec,3) = ncp + 1; % number of nodes per section (nsec)
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,pgrad,1];
                    
                    % knot vector number (ikv)
                    knots(isec_w_NURBS,1) = isec_w_NURBS;
                    % weights number (iw)
                    knots(isec_w_NURBS,2) = isec_w_NURBS;
                    % number of knots (nknots)
                    knots(isec_w_NURBS,3) = nKnot;
                    % initial knot value (iknot)
                    knots(isec_w_NURBS,4) = knotVectors{kvno}(6);
                    % final knot value (jknot)
                    knots(isec_w_NURBS,5) = knotVectors{kvno}(end);
                    % knots
                    knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                    
                    % weights number (iw)
                    wgt(isec_w_NURBS,1) = isec_w_NURBS;
                    % number of control points / weights (nweights)
                    wgt(isec_w_NURBS,2) = ncp;
                    % weights
                    wgt(isec_w_NURBS,3:2+ncp) = coor(elmt(idx(1:end-1)),4);
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
                end
                patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
                hold on
                ii = ii + 1;
            end
        elseif kvno ~= 0 && icp > ecp && ecp == 1
            ii = 1;
            while ii <= nel
                if (ii == icp)  % && (ecp == 1)
                    a = ii;
                    b = 1;
                    ii = ii + ncp-2;
                    idx = [a:nel,1,nel+1];
                    
                    isec = isec + 1;
                    isec_w_NURBS = isec_w_NURBS + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                    sections(isec,3) = ncp + 1; % number of nodes per section (nsec)
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,pgrad,1];
                    
                    
                    % knot vector number (ikv)
                    knots(isec_w_NURBS,1) = isec_w_NURBS;
                    % weights number (iw)
                    knots(isec_w_NURBS,2) = isec_w_NURBS;
                    % number of knots (nknots)
                    knots(isec_w_NURBS,3) = nKnot;
                    % initial knot value (iknot)
                    knots(isec_w_NURBS,4) = knotVectors{kvno}(6);
                    % final knot value (jknot)
                    knots(isec_w_NURBS,5) = knotVectors{kvno}(end);
                    % knots
                    knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                    
                    % weights number (iw)
                    wgt(isec_w_NURBS,1) = isec_w_NURBS;
                    % number of control points / weights (nweights)
                    wgt(isec_w_NURBS,2) = ncp;
                    % weights
                    wgt(isec_w_NURBS,3:2+ncp) = coor(elmt(idx(1:end-1)),4);
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
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
                    isec_w_NURBS = isec_w_NURBS + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                    sections(isec,3) = ncp + 1; % number of nodes per section (nsec)
                    sections(isec,4:4+ncp) = elmt(idx);
                    
                    ord(isec,:) = [isec,pgrad,1];
                    
                   
                    % knot vector number (ikv)
                    knots(isec_w_NURBS,1) = isec_w_NURBS;
                    % weights number (iw)
                    knots(isec_w_NURBS,2) = isec_w_NURBS;
                    % number of knots (nknots)
                    knots(isec_w_NURBS,3) = nKnot;
                    % initial knot value (iknot)
                    knots(isec_w_NURBS,4) = knotVectors{kvno}(6);
                    % final knot value (jknot)
                    knots(isec_w_NURBS,5) = knotVectors{kvno}(end);
                    % knots
                    knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                    
                    % weights number (iw)
                    wgt(isec_w_NURBS,1) = isec_w_NURBS;
                    % number of control points / weights (nweights)
                    wgt(isec_w_NURBS,2) = ncp;
                    % weights
                    wgt(isec_w_NURBS,3:2+ncp) = coor(elmt(idx(1:end-1)),4);
                    
                elseif ii ~= nel
                    a = ii;
                    b = ii+1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
                elseif ii == nel
                    a = nel;
                    b = 1;
                    idx = [a,b,nel+1];
                    
                    isec = isec + 1;
                    sections(isec,1) = isec; % section number (isec)
                    sections(isec,2) = 0; % knot vector number (ikv)
                    sections(isec,3) = 3; % number of nodes per section (nsec)
                    sections(isec,4:6) = elmt(idx);
                    
                    ord(isec,:) = [isec,1,1];
                    
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



