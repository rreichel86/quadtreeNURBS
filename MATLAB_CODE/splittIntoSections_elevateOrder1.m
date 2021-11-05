function [nnode,coor,numsec,maxnsec,sections,ord,knots,wgt,polyElmts,secNQ,secN] = splittIntoSections_elevateOrder1(nnode,coor,numel,connectivity,...
    numKnotVectors,knotVectors,maxnknots,idxControlPoints,seedingPoints_splitt,Quadtree,ep,eq)
% splittIntoSections: Splitt polygonal elements into sections 
%
% INPUT:
% nnnode -------------------- number of coordinates = number of nodes
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
%                              nodes, where the first three entries
%                              iel - element number
%                              ikv - knot vector number
%                              idxLeaf - index of Leaf
%                              which_region - region number
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]
% maxnel --------------------- maximum number of nodes on any element
%
% numKnotVectors ------------- number of knot vectors
% knotVectors ---------------- contains knot vectors and following information
%                              ikv - knot vektor number
%                              degree - NURBS curve degree
%                              icp - index initial control point
%                              ecp - last control point
%                              nkonts - number of knots per knot vector
%
% knotVectors = [ikv, degree, icp, ecp, nknots, knot_1,...,knot_nknots]
% maxnknots ------------------ maximun number of knot on any knot vector
% idxControlPoints ----------- control points indices
% idxControlPoints = [icp, ncp, idx_1,...idx_ncp] 
% Quadtree
% ep/eq ------ elevated p-/q-grad
%
%
%
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
% sections -------------------- sections connectivity matrix as nsec-tupel of 
%                               nodes, where the first three entries
%                               isec - section number
%                               idxLeaf - index of Leaf
%                               ikv - knot vector number
%                               iel - quad element number
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, idxLeaf, ikv, iel, region, nsec, node_1,...,node_nsec]
% 
% ord ------------------------- section polynomial order
% ord = [isec, pgrad, qgrad]
%
% knots = [ikv, iw, nknots, iknot, jknot, knot_1,...,knot_nknots]
% wgt = [iw, nweights, weight_1,...,weigth_nweigths]
%
% polyElmts -------------------- relate sections and polygonal elements
% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
% secNQ -------------------------matrix of all sections from neighbour quad
% secNQ = [idx_sec,numSecNQ,isecNQ1,isecNQ2....isecNQ_numSecNQ]
%
%                     idx_sec --- section number of the unqualified section
%                     numSecNQ --- number of sections from neighbour quad
%                                  per unqualified section
%                     isecNQ --- section number from the neighbour quad 
% 
%
% secN  -------------------------number of the neighbour section
% secN = [idx_sec, isecN]
%                      idx_sec --- section number of the unqualified section
%                      isecN --- number of the neighbour section related to the unqualified section
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
% max number of sections on any polygonal element
maxNumSecPoly = 0;
% list containing the number of control points every NURBS curve segment
% apply function to the 5. and 2. entry of knotVectors cell array. 
% This entries correspond to the degree of the curve and the
% number of knots, respectively.
% ncpoints = (nknot-1) - degree
ncpoints = cellfun(@(x) (x(5)-1)-x(2), knotVectors(1:numKnotVectors));


for ielno = 1:numel
    ikv = connectivity{ielno}(2); % knot vector number 
    nel = connectivity{ielno}(5); % number of node per element
    
    % polygon that have curve edges
    if ikv ~= 0
       ncp = ncpoints(ikv); % number of control points
       numSecPoly = nel - ncp + 2;
       numsec = numsec + nel - ncp + 2; % number of sections
       numsec_w_NURBS = numsec_w_NURBS + 1; % number of sections with 
                                            % boundary defined by a NURBS
       nsec = ncp + 1; % number of nodes per section
       maxncp = max(maxncp, ncp); % max number of control points
       
    % polygon that dont have curve edges   
    else
       numSecPoly = nel;
       numsec = numsec + nel; % number of sections
       nsec = 3; % number of nodes per section
    end
    
    maxNumSecPoly = max(maxNumSecPoly, numSecPoly);
    maxnsec = max(maxnsec,nsec); % max number of nodes on any section 
end

% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
polyElmts = zeros(numel,maxNumSecPoly + 4);
% sections = [isec, idxLeaf, ikv, iel, region, nsec, node_1,...,node_nsec]
sections = zeros(numsec, maxnsec + 6);
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
    idxLeaf = connectivity{ielno}(3); %index of Leaf
    region_nro = connectivity{ielno}(4); % region number 
    nel = connectivity{ielno}(5); % number of nodes per element
    elmt = connectivity{ielno}(6:end); % element connectivity matrix
    ecoor = coor( elmt(1:end), 2:3);
    wg = coor( elmt(1:end), 4);
    
    polyElmts(ielno,1) = ielno;
    polyElmts(ielno,2) = region_nro;
    polyElmts(ielno,3) = nel;
    polyElmts(ielno,end) = idxLeaf;
    
    if kvno ~= 0

        nKnot = knotVectors{kvno}(5);
        pgrad = knotVectors{kvno}(2);
        icp = find( elmt == knotVectors{kvno}(3) );
        ecp = find( elmt == knotVectors{kvno}(4) );

        % number of control points
        ncp = (nKnot - 1) - pgrad;
        
        % update number of sections per polygonal element
        polyElmts(ielno,3) = nel - ncp + 2;

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
    ipsec = 0;
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
            sections(isec,2) = idxLeaf; %index of leaf (idxLeaf)
            sections(isec,3) = 0; % knot vector number (ikv)
            sections(isec,4) = ielno; % element number(iel)
            sections(isec,6) = 3; % number of nodes per section (nsec)
            sections(isec,5) = region_nro; % region number
            sections(isec,7:9) = elmt(idx);

            ord(isec,:) = [isec,1,1];
            
            ipsec = ipsec + 1;
            polyElmts(ielno,3+ipsec) = isec;

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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number(iel)
                sections(isec,6) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:7+ncp) = elmt(idx);

                ord(isec,:) = [isec,pgrad,1];
                
                ipsec = ipsec + 1;
                polyElmts(ielno,3+ipsec) = isec;

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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);


                ord(isec,:) = [isec,1,1];
                
                ipsec = ipsec + 1;
                polyElmts(ielno,3+ipsec) = isec;

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);
                ord(isec,:) = [isec,1,1];
                
            end
            
            ipsec = ipsec + 1;
            polyElmts(ielno,3+ipsec) = isec;
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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number(iel)
                sections(isec,6) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:7+ncp) = elmt(idx);
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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
            
            ipsec = ipsec + 1;
            polyElmts(ielno,3+ipsec) = isec;
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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number(iel)
                sections(isec,6) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:7+ncp) = elmt(idx);

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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
            
            ipsec = ipsec + 1;
            polyElmts(ielno,3+ipsec) = isec;
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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = isec_w_NURBS; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number(iel)
                sections(isec,6) = ncp + 1; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:7+ncp) = elmt(idx);

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
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            elseif ii == nel
                a = nel;
                b = 1;
                idx = [a,b,nel+1];

                isec = isec + 1;
                sections(isec,1) = isec; % section number (isec)
                sections(isec,2) = idxLeaf; %index of leaf
                sections(isec,3) = 0; % knot vector number (ikv)
                sections(isec,4) = ielno; % element number (iel)
                sections(isec,6) = 3; % number of nodes per section (nsec)
                sections(isec,5) = region_nro; % region number
                sections(isec,7:9) = elmt(idx);

                ord(isec,:) = [isec,1,1];

            end
            
            ipsec = ipsec + 1;
            polyElmts(ielno,3+ipsec) = isec;
%                 patch(ecoor(idx,1).', ecoor(idx,2).', color, 'FaceAlpha',.5)
%                 hold on
            ii = ii + 1;
        end
        
    end
end


%% get the neighbours of unqualified sections

%seedingPoints = [isec, isec0, idxLeaf, xcoor, ycoor]

num_seedingPoints = length(seedingPoints_splitt(:,1));


secNQ = cell(num_seedingPoints,1); %array for all section numbers of neighbour quad
%secNQ = {idx_sec,sumSecNQ,sec_NQ1,sec_NQ2, ... ,sec_NQ_sumSecNQ]

secN = []; %array for section numbers of the neighbour section
%secN = [idx_sec,isecN]


for isp = 1: num_seedingPoints
    idx_sec = seedingPoints_splitt(isp,2); %(old) number of unqualified sections 
    ikvo = sections(idx_sec,3);
    
    idxLeaf_sec = seedingPoints_splitt(isp, 3); %idxLeaf of unqualified sections
      
    refLeaf_sec = Quadtree.Node{idxLeaf_sec,1}{2,1}(1:end);
    secNQ{isp}(1,1) = isp;
   
    % search for current Quad neighbours
    % Loop over directions:
    % 1 - West
    % 2 - South
    % 3 - East
    % 4 - North
    idxNQa = [];
    for dir = 1:4 % determine possible neighbour Quad in current direction 
        [exist_NQ, refNQ_sec] = refNeighbour(refLeaf_sec,dir);
        if exist_NQ == 1
            % look for neighbour Quad reference in reference array          
            
            idxNQ = findNeighbour(Quadtree,idxLeaf_sec,refLeaf_sec, refNQ_sec);
            if Quadtree.isleaf(idxNQ)
                idxNQa = [idxNQa;idxNQ];
            else
                idxNQchildren = Quadtree.getchildren(idxNQ);
                    if dir == 1 
                      idxNQc1 = idxNQchildren(3);
                      idxNQc2 = idxNQchildren(4);
                
                    elseif dir == 2  
                      idxNQc1 = idxNQchildren(1);
                      idxNQc2 = idxNQchildren(3); 
                    
                    elseif dir == 3  
                      idxNQc1 = idxNQchildren(1);
                      idxNQc2 = idxNQchildren(2); 
                    elseif dir == 4  
                      idxNQc1 = idxNQchildren(2);
                      idxNQc2 = idxNQchildren(4); 
                                
                    end 
                idxNQa = [idxNQa;idxNQc1;idxNQc2];
            end
        end
    end

    
    % sections = [isec, idxLeaf, ikv, iel,region, nsec, node_1,...,node_nsec]
    % polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]

      %get all section numbers of neighbour quad
      sum_secNQ = 0;
      isecNQ = [];
      
      for j = 1: length(idxNQa)                 
          idxsecNQ = find(idxNQa(j,1) == sections(:,2));
          isecNQ = [isecNQ;idxsecNQ];         
          sum_secNQ = sum_secNQ + length(idxsecNQ);
      end
      secNQ{isp}(1,2) = sum_secNQ;
      secNQ{isp}(1,3:2+length(isecNQ)) = isecNQ';
      
   
      
     % sections = [isec, idxLeaf, ikv, iel,region, nsec, node_1,...,node_nsec]
      %get the neighbour section 
      isecN = 0;
      if ikvo ~= 0 
          sec_NURBS = find(sections(:,3) ~= 0); %get the number of all sections with NURBS curve
          jj = 1;
          while jj <= length(sec_NURBS) && isecN == 0
              idxsecN = sec_NURBS(jj,1);             
              if idxsecN ~= idx_sec
                  if sections(idxsecN,7:end-1) == rot90(sections(idx_sec,7:end-1),2)                 
                      isecN = idxsecN; %section number of the neighbour section 
                  end 
              end
              jj = jj + 1;
          end
          
          
      else %(ikvo == 0)
          
          jj = 1;
          while jj <= length(secNQ{isp}(1,3:end)) && isecN == 0          
              idxsecNQ = secNQ{isp}(1,2+jj);
              if length(sections(idxsecNQ,7:end-1)) == length(sections(idx_sec,7:end-1))                           
                  if sections(idxsecNQ,7:end-2) == rot90(sections(idx_sec,7:end-2),2)            
                      isecN = sections(idxsecNQ,1); %section number of the neighbour section 
                  end 
              end 
              jj = jj + 1;
          end
      end
      secN = [secN;idx_sec,isecN];       
end

%% elevate the p-/q-Grad
[coor,maxnsec,nnode,sections,ord]=elevtGrad(ep,eq,nnode,coor,sections,seedingPoints_splitt,ord,secN,polyElmts,connectivity);

end