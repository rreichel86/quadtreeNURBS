function [coor,nnode,sections,ord] = NodesInsertQ(nnode,coor,sections,ord,polyElmts,connectivity,seedingPoints_splitt,seedingPoints_merge) 


%INPUT:
%
% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 
%
% sections -------------------- sections connectivity matrix as nsec-tupel of 
%                               nodes, where the first six entries
%                               isec - section number
%                               ipoly - polygonal element number
%                               ikv - knot vector number
%                               idxLeaf - index of leaf
%                               region - region number 
%                               nsec - number of nodes per section
% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]
%
% ord = [isec,pgrad,qgrad]
%
% polyElmts -------------------- relate sections and polygonal elements
% polyElmts = [ipoly, region, numSecPoly, sec_1,...,sec_numSecPoly,idxLeaf]
%
%
% connectivity --------------- elements connectivity matrix as nel-tupel of 
%                              nodes, where the first five entries
%                              iel - element number
%                              ikv - knot vector number
%                              idxLeaf - index of Leaf
%                              which_region - region number
%                              nel - number of nodes per element
%
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]
%
%
% seedingPoints = [isec,isec0,idxLeaf,xcoor,ycoor,pgrad,qgrad]
%                         
%                       isec  --  new section number of the (un)qualified section
%                       isec0 --  old section number lof the (un)qualified section 
%                       idxLeaf -- number of Leaf
%                       x/ycoor -- x/y coordinate of the scaling center
%                                  of the unqualified sections
%                       p-/qgrad -- p-/qgrad of section from last calculation


%

%

%
%
%
%OUTPUT: coor, nnode, sections, ord   

% coor = [number, x-coor, y-coor, weight, type, which_region, inside_region]
%
%                              type: 1 -  node
%                                    2 - control point or intersection point
%                                    3 - inserted node
%                              which_region: region number
%                              inside_region: 0 - at the boundary 
%                                             1 - inside 
%                                            -1 - outside 



%% extend sections matrix
% max.qgrad of unqualified sections from last calculation
qgrad_splitt = max(seedingPoints_splitt(:,8));

if isempty(seedingPoints_merge) == 0
    % max.qgrad of unqualified sections from last calculation
    qgrad_merge = max(seedingPoints_merge(:,8));
else
    qgrad_merge = 1;
end

% max.pgrad of all sections
maxpgrad = max(ord(:,2));
%number of sections
numsec = length(sections(:,1));

if qgrad_splitt + 1 <= qgrad_merge 
    sections = [sections,zeros(numsec,(qgrad_merge-1)*(maxpgrad+1))];
else
    sections = [sections,zeros(numsec,qgrad_splitt*(maxpgrad+1))];
end

          
%% insert nodes in Q-direction for unqualified sections

% sections = [isec, ipoly, idxLeaf, ikv, region, nsec, node_1,...,node_nsec]

numSeedingPoints_splitt = length(seedingPoints_splitt(:,1));

ElmtUQsec = [];    %element which contain unqualified sections
for isp = 1: numSeedingPoints_splitt
    isec0 = seedingPoints_splitt(isp,2);
    qgrad = seedingPoints_splitt(isp,8);
    iel_splitt = sections(isec0,2);
    ElmtUQsec = [ElmtUQsec;iel_splitt,qgrad];
end
[~,ia] = unique(ElmtUQsec(:,1));
ElmtUQsec = ElmtUQsec(ia,:);
    
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]
[numElmtUQsec,~] = size(ElmtUQsec); %number of poly element
for i = 1: numElmtUQsec
    iel = ElmtUQsec(i); %element number
    qgrad = ElmtUQsec(i,2); %qgrad from last calculation
    elmt = connectivity{iel}(1,6:end); % element connectivity matrix
    kvno = connectivity{iel}(2); %knot vector for element
    numSecPoly = polyElmts(iel,3);
    secElmts = polyElmts(iel,4:3 + numSecPoly);
    eq = qgrad + 1; %elevated qgrad
    ninode = eq - 1; %number of inserted nodes
    
    if ninode == 0
        break
    end
       
    
    if kvno == 0 
       
        
        %insert nodes in shared edges between sections per element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = []; %vertices of this element
        for inelmt = 1:length(elmt)-1
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end
        
     
               
        %insert nodes per section within element  
             
        for ise = 1:numSecPoly
            idxse = secElmts(ise);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = []; %nodes of each section except vertices of element
            %insert nodes 
            for indsec = 1:nsec-1                        
                inode = sections(idxse,6+indsec);
                if ismember(inode,nodeElmt(:,1)) == 0 % vertices of element are not including
                    xcoor1 = coor(inode,2);
                    ycoor1 = coor(inode,3);
                    coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                    %update coor matrix
                    coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                    coor(nnode+1:nnode+ninode,4) = 1;  %weight
                    coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                    coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                    coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
                    coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                    nodeSec0 =[inode,nnode+1:nnode + ninode];
                    nodeSec = [nodeSec;nodeSec0];
                    nnode = nnode + ninode; 
                end
            end
                      
            
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center
            
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
        
    else %kvno ~= 0
        
        %insert nodes in shared edges between sections per element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = [];
        for inelmt = 1:length(elmt)-1
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end    
                
  
        %insert nodes per section within element       
        
        for ise = 1:numSecPoly
            idxse = secElmts(ise); %section number within element
            ikv = sections(idxse,4);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = [];
            %insert nodes 
            if ikv == 0                
                for indsec = 1:nsec-1                        
                    inode = sections(idxse,6+indsec);
                    if ismember(inode,nodeElmt(:,1)) == 0
                        xcoor1 = coor(inode,2);
                        ycoor1 = coor(inode,3);
                        coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                        %update coor matrix
                        coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                        coor(nnode+1:nnode+ninode,4) = 1;  %weight
                        coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                        coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                        coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
                        coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                        nodeSec0 =[inode,nnode+1:nnode + ninode];
                        nodeSec = [nodeSec;nodeSec0];
                        nnode = nnode + ninode; 
                    end
                end 
            end
            
         
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center                
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
    end
end

%% insert nodes in Q-direction for qualified sections

if isempty(seedingPoints_merge) == 1
    return
end

numSeedingPoints_merge = length(seedingPoints_merge(:,1));

ElmtQsec = [];    %element which contain qualified sections
for isp = 1: numSeedingPoints_merge
    isec0 = seedingPoints_merge(isp,2);
    qgrad = seedingPoints_merge(isp,8);
    iel_merge = sections(isec0,2);
    if ismember(iel_merge,ElmtUQsec(:,1)) == 0
        ElmtQsec = [ElmtQsec;iel_merge,qgrad];
    end
end

[~,ia] = unique(ElmtQsec(:,1));
ElmtQsec = ElmtQsec(ia,:);

if isempty(ElmtQsec) == 1
    return
end

    
% connectivity = [iel, ikv, idxLeaf, which_region, nel, node_1,...,node_nel, scaling_center]
[numElmtQsec,~] = size(ElmtQsec); %number of poly element
for i = 1: numElmtQsec
    iel = ElmtQsec(i,1); %element number
    qgrad = ElmtQsec(i,2); % qgrad from last calculation
    elmt = connectivity{iel}(1,6:end); % vertices of element and scaling_center 
    kvno = connectivity{iel}(2); %knot vector for element
    numSecPoly = polyElmts(iel,3);
    secElmts = polyElmts(iel,4:3 + numSecPoly);
    eq = qgrad; %elevated qgrad 
    ninode = eq - 1; %number of inserted nodes
    
    if ninode == 0
        break
    end
    
    if kvno == 0 
       
        
        %insert nodes in shared edges between sections per element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = []; %vertices of this element
        for inelmt = 1:length(elmt)-1
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end
        
     
               
        %insert nodes per section within element       
        
        for ise = 1:numSecPoly
            idxse = secElmts(ise);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = []; %nodes of each section except vertices of element
            %insert nodes 
            for indsec = 1:nsec-1                        
                inode = sections(idxse,6+indsec);
                if ismember(inode,nodeElmt(:,1)) == 0 % vertices of element are not including
                    xcoor1 = coor(inode,2);
                    ycoor1 = coor(inode,3);
                    coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                    %update coor matrix
                    coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                    coor(nnode+1:nnode+ninode,4) = 1;  %weight
                    coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                    coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                    coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
                    coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                    nodeSec0 =[inode,nnode+1:nnode + ninode];
                    nodeSec = [nodeSec;nodeSec0];
                    nnode = nnode + ninode; 
                end
            end
                      
            
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center
            
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
        
    else %kvno ~= 0
        
        %insert nodes in shared edges between sections per element
        xcoor2 = coor(elmt(end),2); %x_coordinate of sc
        ycoor2 = coor(elmt(end),3); %y_coordinate of sc  
        nodeElmt = [];
        for inelmt = 1:length(elmt)-1
            inode = elmt(inelmt);
            xcoor1 = coor(inode,2);
            ycoor1 = coor(inode,3);
            coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
            %update coor matrix
            coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
            coor(nnode+1:nnode+ninode,4) = 1;  %weight
            coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
            coor(nnode+1:nnode+ninode,6) = 0;  %which-region
            coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
            coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
            nodeElmt0 =[inode,nnode+1:nnode + ninode];
            nodeElmt = [nodeElmt;nodeElmt0];
            nnode = nnode + ninode;
        end    
                
  
        %insert nodes per section within element       
        
        for ise = 1:numSecPoly
            idxse = secElmts(ise); %section number within element
            ikv = sections(idxse,4);
            nsec = sections(idxse,6);
            pgrad = ord(idxse,2);
            
            nodeSec = [];
            %insert nodes 
            if ikv == 0                
                for indsec = 1:nsec-1                        
                    inode = sections(idxse,6+indsec);
                    if ismember(inode,nodeElmt(:,1)) == 0
                        xcoor1 = coor(inode,2);
                        ycoor1 = coor(inode,3);
                        coor_inodes_q = getLocation(xcoor1,xcoor2,ycoor1,ycoor2,eq);
                        %update coor matrix
                        coor(nnode+1:nnode+ninode,1) = nnode+1: nnode + ninode;
                        coor(nnode+1:nnode+ninode,4) = 1;  %weight
                        coor(nnode+1:nnode+ninode,5) = 3;  %typ(inserted nodes)
                        coor(nnode+1:nnode+ninode,6) = 0;  %which-region
                        coor(nnode+1:nnode+ninode,7) = -1; %inside-region  
                        coor(nnode+1:nnode+ninode,2:3) =[coor_inodes_q];
                        nodeSec0 =[inode,nnode+1:nnode + ninode];
                        nodeSec = [nodeSec;nodeSec0];
                        nnode = nnode + ninode; 
                    end
                end 
            end
            
         
            
            %divide inserted nodes
            sc = sections(idxse,6+nsec); %number of scaling center                
            for indsec =1: nsec-1
                inode = sections(idxse,6+indsec);
                for eta = 1:ninode                   
                    if ismember(inode,nodeElmt(:,1)) == 1
                        idx = find(inode == nodeElmt(:,1));
                        idxnode = nodeElmt(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    else 
                        idx = find(inode == nodeSec(:,1));
                        idxnode = nodeSec(idx,eta+1);
                        sections(idxse,6+indsec+(pgrad+1)*eta) = idxnode;
                    end
                end
            end
            
            sections(idxse,6) = (pgrad+1)*eq + 1; %nsec
            sections(idxse,6+(pgrad+1)*eq + 1) = sc; %scaling center
            
            ord(idxse,3) = eq;
        end
    end
end 
    
%% adjust the size of sections matrix
maxnsec = max(sections(:,6));
sections(:,6+maxnsec+1:end) = [];            

end
                                
    