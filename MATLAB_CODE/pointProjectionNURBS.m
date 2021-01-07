function [Q,v] = pointProjectionNURBS(P,u0,NURBS)


u = zeros(11,1);

eps1 = 1e-12;
eps2 = 1e-12;

a = NURBS.knots(1);
b = NURBS.knots(end);

u(1) = u0;

for itr = 1 : 10
    
    [R,R_xi,R2_xi] = shape_func_NURBS_1d_2(u(itr),NURBS.degree,NURBS.knots,NURBS.weights);
    
    C = NURBS.controlPoints*R;
    dCdu = NURBS.controlPoints*R_xi;
    dC2du = NURBS.controlPoints*R2_xi;
    
    f = dCdu'*(C - P);
    dfdu = dC2du'*(C - P) + norm(dCdu)^2;
    
    % point coincidence
    cond1 = norm(C - P) < eps1;
    % zero cosine
    cond2 = norm(f)/norm(dCdu)/norm(C-P) < eps2;
    
    if cond1 || cond2
        Q = C;
        v = u(itr);
        break
    end
    
    u(itr+1) =  u(itr) - f/dfdu;
    
    % ensure that the parameter stays within the range
    if u(itr+1) < a
        u(itr+1) = a;
    end
    
    if u(itr+1) > b
        u(itr+1) = b;
    end
    
    % parameter does not change significantly
    du = u(itr+1) - u(itr);
    cond4 = norm( du*dCdu ) < eps1;
    
    if cond4
        Q = C;
        v = u(itr);
        break
    end
    
end



end