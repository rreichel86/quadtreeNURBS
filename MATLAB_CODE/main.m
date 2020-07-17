clear all;
close all;
clc;

% Select example
example_nro = 2;
%      Nr.     |      Description        |    Section     
%---------------------------------------------------------
%       1      |       Moby-Dick         |      3         
%       2      |     Circumference       |     4.1
%       3      |    Double circumf.      |     4.1
%       4      |      Flat shape         |     4.1
%       5      |  Limit double circum.   |     4.2
%       6      |    Degenerated case     |     4.2



% Initialization
% ==============
% Control points input should be a matrix with dimensions
% (coordinates,nPoints), and should be as follows: first row corresponds to
% the x coordinate, second to the y coordinate, . Each control point
% consequently lies in its own column.

% Obtains the selected NURBS definition
[degree,knots,controlPoints,weights,ax,Boundary]=NURBS_parameters(example_nro);

% CalculateNURBS computes the curve
NURBS=CalculateNURBS(degree,knots,controlPoints,weights);

% First plotting
figure(1);
hold on;
axis(ax);
plot(NURBS(1,:),NURBS(2,:),'r','LineWidth',3);
plot(controlPoints(1, :), controlPoints(2, :), 'ro','LineWidth',3);
plot(controlPoints(1, :), controlPoints(2, :), '--','LineWidth',0.5);
box on



%Count of control points in the root
nPoints = checkQuad( Boundary(:,2)',  Boundary(:,4)', Boundary(:,2)',...
                    controlPoints);

% Split the root if there is more than one point in the root
if nPoints > 1
    %Setting tree
    Quadtree = tree('root');
    data=Quadtree.Node{1,1};
    data={data;[]};
    Quadtree = Quadtree.set(1, data);
    
    l=1; % l: 1 if quad NW, 2 SW, 3 NE, 4 SE. Initialize at 1
    k=0; % k: level of decomposition, equal to 0 at the root
    pos_aux=[]; Q_aux=[0]; %auxiliar arrays
        
    [Quadtree] = decompose(Quadtree,Boundary,NURBS ,controlPoints, knots,...
                            weights, degree, l, k,pos_aux, Q_aux,Boundary);
end

%hold off
                        
[Quadtree]=Star_Shape(Quadtree,controlPoints,NURBS, knots, weights,...
                        degree, Boundary);


[Quadtree]=Balance_Quadtree(Quadtree,NURBS,controlPoints, knots, ...
                            weights, degree,Boundary);

% close all the figures and plots the NURBS contained in each leaf separately
%plot_leaf(Quadtree,ax)
[coor,connectivity,nnode,maxnel,numel,kv_element,kv_num,maxnk]=extract_leaf(Quadtree);

%% Plot polygonal elements 
 for ielno = 1:numel
     nel = connectivity{ielno}(2);
     elmt = connectivity{ielno}(3:2+nel);
     ex = coor( elmt, 2);
     ey = coor( elmt, 3);
     
     patch(ex,ey, 'red','FaceAlpha',.5)
 end 


%% Plot polygonal elements with curve edges
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
        
        iknot_coor = ecoor(iknot,:);
        eknot_coor = ecoor(eknot,:);
        
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
%          continue
%          for ii = 1:nel
%              if ii ~= nel 
%                  a = ii;
%                  b = ii+1;
%              else
%                  a = nel;
%                  b = 1;
%              end    
%              plot(ecoor([a,b],1).', ecoor([a,b],2).', 'c-')
%              
%              hold on 
%          end   
         
         patch(ecoor(:,1),ecoor(:,2), 'green','FaceAlpha',.5)
         
      

     % plot polygon that dont have curve edges
     elseif kvno ~= 0 && iknot < eknot && iknot == 1
%          continue
         ii = 1;
         while ii <= nel 
            if (ii == eknot) && (iknot == 1)
                 a = ii;
                 b = 1;
                 ii = ii + pgrad - 1;
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
%          continue
         ii = 1;
         while ii <= nel 
             if ii == iknot
                 a = ii;
                 b = ii + pgrad;
                 ii = ii + pgrad - 1;
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
%          continue
         ii = 1;
         while ii <= nel 
             if (ii == iknot) && (eknot == 1)
                 a = ii;
                 b = 1;
                 ii = ii + pgrad - 1;
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
%          continue
         ii = 1;
         while ii <= nel 
             if ii == eknot
                 a = ii;
                  b = ii + pgrad;
                 ii = ii + pgrad - 1;
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



