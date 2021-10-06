function [nnode,coor,numsec,maxnsec,sections,ord,knots,wgt] = splittIntoSections_2(nnode,coor,numel,connectivity,...
    numKnotVectors,knotVectors,idxControlPoints)
% splittIntoSections: Splitt polygonal elements into section
%
% INPUT:
% nnode ---------------------- number of coordinates = number of nodes
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
% numel ---------------------- number of elements
% connectivity --------------- elements connectivity matrix as nel-tupel of 
%                              nodes, where the first six entries
%                              iel - element number
%                              ikv - knot vector number 
%                              which_region - region number
%                              inode - index initial node 
%                              jnode - index last node
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, which_region, inode, jnode, nel, node_1,...,node_nel, scaling_center]
% maxnel --------------------- maximum number of nodes on any element
%
% numKnotVectors ------------- number of knot vectors
% knotVectors ---------------- contains knot vectors and following information
%                              ikv - knot vektor number
%                              degree - NURBS curve degree
%                              iknot - initial knot value
%                              jknot - end knot value
%                              nkonts - number of knots per knot vector
%
% knotVectors = [ikv, degree, iknot, jknot, nknots, knot_1,...,knot_nknots]
% idxControlPoints ----------- control points indices
% idxControlPoints = [icp, ncp, idx_1,...idx_ncp] 
%
% OUTPUT:
% nnode ---------------------- number of coordinates = number of nodes
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
% ord = [isec, pgrad, qgrad, pos, Npos]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]
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
% number of sections which boundary is defined by a NURBS
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
    icp = connectivity{ielno}(1,4);
    ecp = connectivity{ielno}(1,5);
    
    % polygon that have curve edges
    if ikv ~= 0
        ncp = ncpoints(ikv); % number of control points
        numsec = numsec + nel; % number of sections
        if icp == ecp 
            numsec_w_NURBS = numsec_w_NURBS + 3; % number of sections with
        else
            numsec_w_NURBS = numsec_w_NURBS + 2;
        end     
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
        
        % control points indices
        idxCtrlP = idxControlPoints{kvno}(3:end);
        
        NURBS_segment = CalculateNURBS_2(pgrad,iKnot,jKnot,knotVectors{kvno}(6:end),...
            coor(idxCtrlP,2:3)',coor(idxCtrlP,4)');
        
        OP = orientation( ecoor(end,:), NURBS_segment(1,1:2) ,NURBS_segment(end,1:2) );
        
    end
    
    % polygon that dont have curve edges
    if kvno == 0
        for ii = 1:nel
            a = ii;
            if ii ~= nel
                b = ii+1;
            else
                b = 1;
            end
            
            idx = [a,b,nel+1];
            
            isec = isec + 1;
            sections(isec,1:4) = [isec,0,region_nro,3];
            sections(isec,5:7) = elmt(idx);
            ord(isec,:) = [isec,1,1,0,0];
        end
    % polygon that have curve edges
    elseif kvno ~= 0  % && icp ~= ecp   % && icp ~= 1
        ii = 1;
        while ii <= nel
            if ii == icp  && icp ~= ecp% section with a NURBS segment 
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                
                sections(isec,1:4) = [isec,isec_w_NURBS,region_nro,ncp + 1];
                if OP == -1
                    sections(isec,5:4+ncp) = fliplr(idxCtrlP);
                else
                    sections(isec,5:4+ncp) = idxCtrlP;
                end
                sections(isec,5+ncp) = elmt(end);
                
                ord(isec,:) = [isec,pgrad,1,0,0];
                
                % knots
                knots(isec_w_NURBS,1:3) = [isec_w_NURBS,isec_w_NURBS,nKnot];
                if OP == -1
                    knots(isec_w_NURBS,4:5) = [1-jKnot,1-iKnot];
                else
                    knots(isec_w_NURBS,4:5) = [iKnot,jKnot];
                end
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights
                wgt(isec_w_NURBS,1:2) = [isec_w_NURBS,ncp];
                if OP == -1
                    wgt(isec_w_NURBS,3:2+ncp) = coor(fliplr(idxCtrlP),4);
                else
                    wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                end
                
            elseif ( ii == (icp - 1) && icp ~= 1 ) ||  ( ii == nel && icp == 1 ) % section adjacent to a NURBS segment 
                a = ii;
                b = icp;
               
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1:4) = [isec,isec_w_NURBS,region_nro,ncp + 2];
                
                sections(isec,5) = elmt(a);
                sections(isec,6:5+ncp) = idxCtrlP;
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,1,0];
                if norm(coor(elmt(b),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knots
                knots(isec_w_NURBS,1:5) = [isec_w_NURBS,isec_w_NURBS,nKnot,iKnot,jKnot];
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights
                wgt(isec_w_NURBS,1:2) = [isec_w_NURBS,ncp];
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
            elseif ii == ecp % section adjacent to a NURBS segment 
                a = ii;
                if ii == nel
                    b = 1;
                else
                    b = ii + 1;
                end
                
                isec = isec + 1;
                isec_w_NURBS = isec_w_NURBS + 1;
                sections(isec,1:4) = [isec,isec_w_NURBS,region_nro,ncp + 2];
                
                sections(isec,5:4+ncp) = idxCtrlP;
                sections(isec,5+ncp) = elmt(b);
                sections(isec,6+ncp) = elmt(end);
                
                ord(isec,:) = [isec,1,1,-1,0];
                if norm(coor(elmt(a),2:3) - NURBS_segment(1,1:2)) < 1e-10
                    ord(isec,5) = -1;
                else
                    ord(isec,5) = 1;
                end
                
                % knots
                knots(isec_w_NURBS,1:5) = [isec_w_NURBS,isec_w_NURBS,nKnot,iKnot,jKnot];
                knots(isec_w_NURBS,6:5+nKnot) = knotVectors{kvno}(6:end);
                
                % weights
                wgt(isec_w_NURBS,1:2) = [isec_w_NURBS,ncp];
                wgt(isec_w_NURBS,3:2+ncp) = coor(idxCtrlP,4);
                
            elseif ii <= nel
                a = ii;
                if ii == nel
                    b = 1;
                else
                    b = ii+1;
                end
                
                idx = [a,b,nel+1];
                
                isec = isec + 1;
                sections(isec,1:4) = [isec,0,region_nro,3];
                sections(isec,5:7) = elmt(idx);
                ord(isec,:) = [isec,1,1,0,0];
            end
            ii = ii + 1;
        end
    else 
        kvno
        icp 
        ecp 
    end
end

for isec = 1:numsec
    nsec = sections(isec,4);
    elmt = sections(isec,5:4+nsec)';
    coor(elmt,1) = 1;
end 

idx = 0;
for i = 1:nnode
   
    if coor(i,1) == 1
        idx = idx + 1;
        coor(i,1) = idx;
    else
        coor(i,1) = -99;
    end     
    
end    

% update element connectivity 
for isec = 1:numsec
    nsec = sections(isec,4);
    new_elmt = coor(sections(isec,5:4+nsec), 1)';
    sections(isec,5:4+nsec) = new_elmt;
end 

nnode = idx;
coor = coor(coor(:,1) ~= -99,:);



end