function [nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = nurbs_quadtree_mesh(k_min,degree,knots,controlPoints,weights,Boundary)

% OUTPUT:
% nnode -------------------- number of coordinates = number of nodes
% coor ----------------------- nodes coordinates and weights 
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary
%                                             1 - inside 
%                                            -1 - outside
%
% numsec ---------------------- number of sections
% maxnsec --------------------- maximum number of nodes on any section
% sections -------------------- sectionsconnectivity matrix as nsec-tupel of 
%                               nodes, where the first three entries
%                               isec - section number
%                               ikv - knot vector number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ikv, region, nsec, node_1,...,node_nsec]
% 
% ord ------------------------- section polynomial order
% ord = [isec, pgrad, qgrad]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]


%% Quadtree decomposition
[Quadtree] = nurbs_brep_quadtree(k_min,degree,knots,controlPoints,weights,Boundary);

%% Extract polygonal elements 
[nnode,coor,numel,connectivity,maxnel,...
 numKnotVectors,knotVectors,maxnknots,idxControlPoints] = extractElements(Quadtree);

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
    nel = connectivity{ielno}(4);
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

% sections = [isec, ikv, region, nsec, node_1,...,node_nsec]
sections = zeros(numsec, maxnsec + 4);
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
knots = zeros(numsec_w_NURBS, maxnknots + 5);
% ord = [isec, pgrad, qgrad]
ord = zeros(numsec, 3);
% weights = [iw, nweights, weight_1,...,weigth_nweigths]
wgt = zeros(numsec_w_NURBS, maxncp + 2);

isec = 0;
isec_w_NURBS = 0;
for ielno = 1:numel
    kvno = connectivity{ielno}(2); % knot vector number
    region_nro = connectivity{ielno}(3); % region number 
    nel = connectivity{ielno}(4); % number of nodes per element
    elmt = connectivity{ielno}(5:end); % element connectivity matrix
    ecoor = coor( elmt(1:end), 2:3);
    wg = coor( elmt(1:end), 4);

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
            OP = 1;
        else
            OP = -1;
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
            sections(isec,4) = 3; % number of nodes per section (nsec)
            sections(isec,3) = region_nro; % region number
            sections(isec,5:7) = elmt(idx); 

            ord(isec,:) = [isec,1,1];

%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
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
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:5+ncp) = elmt(idx);

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
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = 0; % knot vector number (ikv)
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
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
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:5+ncp) = elmt(idx);

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
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = 0; % knot vector number (ikv)
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
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
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:5+ncp) = elmt(idx);

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
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = 0; % knot vector number (ikv)
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
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
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:5+ncp) = elmt(idx);

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
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = 0; % knot vector number (ikv)
                sections(isec,4) = 3; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:7) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
            ii = ii + 1;
        end
    end
end



end
