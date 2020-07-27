function [R,dR]=Nurbs_shape_function(u,p,U,w);

%Input

%u:parameterical point where Basis functions and it's derivatives have
%to find.
%p:degree of the basis functions
%cp:control points of NURBS curve
%U:Knot vector of NURBS curve
%w:Weights of control points 

%Output

%N:NURBS basis function values
%N_x:NURBS derivative function


n=length(U)-p-1;

s = findspan(n,p,u,U);
[N_basis,dN_basis] = basis_DerFun(s,u,p,U)
s=s+1;
w=[w(1,s-p:s)];
sum_N_basis_w=sum(N_basis.*w);
sum_dN_basis_w=sum(dN_basis.*w);

for i=1:p+1
    R(1,i)=(N_basis(1,i)*w(1,i))/sum_N_basis_w;
    dR(1,i)=((dN_basis(1,i)*w(1,i))/sum_N_basis_w)-(N_basis(1,i)*w(1,i)*sum_dN_basis_w/(sum_N_basis_w)^2);
end
end
    


