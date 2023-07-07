function plotNURBSCurve(NURBS)

% compute point of the NURBS curve
NURBS_pts = CalculateNURBS(NURBS);

% Plot control polygon
plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'b-.','LineWidth',1);
% Plot control points
plot(NURBS.controlPoints(1, :), NURBS.controlPoints(2, :),'o','Color','r','MarkerFaceColor','r','MarkerSize',6);
% Plot curve
plot(NURBS_pts(:,1),NURBS_pts(:,2),'m','LineWidth',2.0);



end