function [nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = splittIntoSections_2(nnode,coor,numel,connectivity,...
    numKnotVectors,knotVectors,idxControlPoints,indices)
% splittIntoSections: Splitt polygonal elements into section
%
% INPUT:
% coor
% numel
% connectivity
% numKnotVectors
% knotVectors
% maxnknots
% idxCotrolPoints
%
% OUTPUT:
% coor
% numsec
% numsec_w_NURBS
% maxnsec
% maxncp
% sections
% knots
% ord
% wgt
%
% -----------------------------------------------------------------------------

% Compute:
% numsec         --- number of sections
% numsec_w_NURBS --- number of sections with boundary defined by a NURBS
% maxnsec        --- max number of node on any section
% maxncp         --- max number of control points on any section
%                    with boundary defined by a NURBS
% ncpoints       --- list containing the number of control points every
%                    NURBS curve segment

% number of sections
numsec = 0;
% number of sections with boundary defined by a NURBS
numsec_w_NURBS = 0;
% max number of nodes on any section
maxnsec = 0;
% max number of control points
maxncp = 0;
% list containing the number of control points every NURBS curve segment
% apply function to the 5. and 2. entry of knotVectors cell array.
% This entries correspond to the degree of the curve and the
% number of knots, respectively.
% ncpoints = (nknot-1) - degree
ncpoints = cellfun(@(x) (x(5)-1)-x(2), knotVectors(1:numKnotVectors));
nknots = cellfun(@(x) x(5), knotVectors(1:numKnotVectors));
maxnknots = max(nknots);

for ielno = 1:numel
    ikv = connectivity{ielno}(2); % knot vector number
    nel = connectivity{ielno}(6); % number of node per element
    
    % polygon that have curve edges
    if ikv ~= 0
        ncp = ncpoints(ikv); % number of control points
        numsec = numsec + nel; % number of sections
        numsec_w_NURBS = numsec_w_NURBS + 3; % number of sections with
        % boundary defined by a NURBS
        nsec = ncp + 1; % number of nodes per section
        maxncp = max(maxncp, ncp); % max number of control points
        
        % polygon that dont have curve edges
    else
        numsec = numsec + nel; % number of sections
        nsec = 3; % number of nodes per section
    end
    
    maxnsec = max(maxnsec,nsec); % max number of nodes on any section
end

% sections = [isec, ikv, region, nsec, node_1,...,node_nsec]
sections = zeros(numsec, maxnsec + 4);
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
knots = zeros(numsec_w_NURBS, maxnknots + 5);
% ord = [isec, pgrad, qgrad]
ord = zeros(numsec, 5);
% weights = [iw, nweights, weight_1,...,weigth_nweigths]
wgt = zeros(numsec_w_NURBS, maxncp + 2);

isec = 0;
isec_w_NURBS = 0;

for ielno = 1:numel
    kvno = connectivity{ielno}(2); % knot vector number
    region_nro = connectivity{ielno}(3); % region number
    nel = connectivity{ielno}(6); % number of nodes per element
    elmt = connectivity{ielno}(7:end); % element connectivity matrix
    ecoor = coor( elmt(1:end), 2:3);
    
    if kvno ~= 0
        
        nKnot = knotVectors{kvno}(5);
        pgrad = knotVectors{kvno}(2);
        iKnot = knotVectors{kvno}(3);
        jKnot = knotVectors{kvno}(4);
        
        icp = connectivity{ielno}(1,4);
        ecp = connectivity{ielno}(1,5);
        
        % number of control points
        ncp = (nKnot - 1) - pgrad;
        
        % scaling center coordinates
        sc_coor = ecoor(end,:);
        
        % control points indices
        idxCtrlP = idxControlPoints{kvno}(3:end);
        
        NURBS_segment = CalculateNURBS_2(pgrad,iKnot,jKnot,knotVectors{kvno}(6:end),...
            coor(idxCtrlP,2:3)',coor(idxCtrlP,4)');
        
        nSegPts = size(NURBS_segment,1);
        
        poly = zeros(nSegPts+1,2);
        
        poly(1:nSegPts,:) = NURBS_segment(:,1:2);
        poly(nSegPts+1,:) = sc_coor;
        
        % determine signed polygon area
        % area =  1 for CCW
        % area = -1 for CW
        area = poly(nSegPts+1,1) * poly(1,2) - poly(1,1) * poly(nSegPts+1,2);
        
        for i = 1: nSegPts
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
            
            ord(isec,:) = [isec,1,1,0,0];
            
        end
        
        
        
        % plot polygon that have curve edges
    elseif kvno ~= 0 && icp < ecp && icp ~= 1
        ii = 1;
        while ii <= nel
            if ii == icp
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
%                 if OP == -1
%                     sections(isec,5:4+ncp) = fliplr(idxCtrlP);
%                 else
                    sections(isec,5:4+ncp) = idxCtrlP;
%                 end
                sections(isec,5+ncp) = elmt(end);
                
                ord(isec,:) = [isec,pgrad,1,0,0];
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                
                
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                
                
                
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
%                 if OP == -1
%                     wgt(isec_w_NURBS,3:2+ncp) = coor(fliplr(idxCtrlP),4);
%                 else
                    wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
%                 end
                
            elseif ii == icp - 1
                a = ii;
                b = ii+1;
                %idx = [a,b,nel+1];
                
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 2; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5) = elmt(a);
                sections(isec,6:5+ncp) = idxCtrlP;
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,1,0];
                
                if norm(coor(elmt(b),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
            elseif ii == ecp
                a = ii;
                if ii == nel
                    b = 1;
                else
                    b = ii + 1;
                end
                %idx = [a,b,nel+1];
                
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 2; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:4+ncp) = idxCtrlP;
                sections(isec,5+ncp) = elmt(b);
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,-1,0];
                
                if norm(coor(elmt(a),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
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
                
                ord(isec,:) = [isec,1,1,0,0];
                
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
                
                ord(isec,:) = [isec,1,1,0,0];
                
            end
            ii = ii + 1;
        end
    elseif kvno ~= 0 && icp > ecp && ecp == 1
        ii = 1;
        while ii <= nel
            if (ii == icp)  % && (ecp == 1)
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
%                 if OP == -1
%                     sections(isec,5:4+ncp) = fliplr(idxCtrlP);
%                 else
                    sections(isec,5:4+ncp) = idxCtrlP;
%                 end
                sections(isec,5+ncp) = elmt(end);
                
                ord(isec,:) = [isec,pgrad,1,0,0];
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                
                
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
%                 if OP == -1
%                     wgt(isec_w_NURBS,3:2+ncp) = coor(fliplr(idxCtrlP),4);
%                 else
                    wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
%                 end
                
            elseif ii == icp - 1
                a = ii;
                b = ii+1;
                %idx = [a,b,nel+1];
                
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 2; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5) = elmt(a);
                sections(isec,6:5+ncp) = idxCtrlP;
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,1,0];
                
                if norm(coor(elmt(b),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
            elseif  ii == ecp
                a = ii;
                b = ii+1;
                %idx = [a,b,nel+1];
                
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ncp + 2; % number of nodes per section (nsec)
                sections(isec,3) = region_nro; % region number
                sections(isec,5:4+ncp) = idxCtrlP;
                sections(isec,5+ncp) = elmt(b);
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,-1,0];
                
                if norm(coor(elmt(a),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knot vector number (ikv)
                knots(isec_w_NURBS,1) = isec_w_NURBS;
                % weights number (iw)
                knots(isec_w_NURBS,2) = isec_w_NURBS;
                % number of knots (nknots)
                knots(isec_w_NURBS,3) = nKnot;
                % initial knot value (iknot)
                knots(isec_w_NURBS,4) = iKnot;
                % final knot value (jknot)
                knots(isec_w_NURBS,5) = jKnot;
                % knots
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights number (iw)
                wgt(isec_w_NURBS,1) = isec_w_NURBS;
                % number of control points / weights (nweights)
                wgt(isec_w_NURBS,2) = ncp;
                % weights
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
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
                
                ord(isec,:) = [isec,1,1,0,0];
                
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
                
                ord(isec,:) = [isec,1,1,0,0];
                
            end
            ii = ii + 1;
        end
    end
end

end