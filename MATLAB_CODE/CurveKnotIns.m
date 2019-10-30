function [knots_new, Q, weights] = CurveKnotIns(degree, controlPoints, knots, weights, newKnot)
% Adapted form the NURBS book
p=degree;
np = length(controlPoints);
s=0;
r=1;

k=find(knots>=newKnot(1),1)-1;

Pw=[controlPoints(1,:) .* weights ;controlPoints(2,:) .* weights; weights];

mp= np+degree+1;
nq= np+ r;

for i = 1:k 
    knots_new(i) = knots(i);
end

for i = 1:r
knots_new(k+i)= newKnot;
end

for i = k+1:mp 
    knots_new(i+r) = knots(i);

end

for i = 1:k-p 
    Qw(1,i) = Pw(1,i);
    Qw(2,i) = Pw(2,i);
    Qw(3,i) = Pw(3,i);

end

for i = k-s:np
    Qw(1,i+r) = Pw(1,i);
    Qw(2,i+r) = Pw(2,i);
    Qw(3,i+r) = Pw(3,i);
end

for i = 1:p-s+1 
    Rw(1,i) = Pw(1,k-p+i-1);
    Rw(2,i) = Pw(2,k-p+i-1);
    Rw(3,i) = Pw(3,k-p+i-1);

end

for j = 1:r
    L= k-p+j;
    

    for i = 1:p-j-s+1
        alpha = (newKnot- knots(L+i-1))/(knots(i+k) - knots(L+i-1));
        Rw(1,i) = alpha*Rw(1,i+1) + (1.0-alpha)*Rw(1,i); 
        Rw(2,i) = alpha*Rw(2,i+1) + (1.0-alpha)*Rw(2,i); 
        Rw(3,i) = alpha*Rw(3,i+1) + (1.0-alpha)*Rw(3,i); 
    end
    
    Qw(1,L) = Rw(1,1);
    Qw(2,L) = Rw(2,1);
    Qw(3,L) = Rw(3,1);
    Qw(1,k+r-j-s) = Rw(1,p-j-s+1); 
    Qw(2,k+r-j-s) = Rw(2,p-j-s+1);
    Qw(3,k+r-j-s) = Rw(3,p-j-s+1);
end

for i = (L+1):(k-s-1) 
    Qw(1,i) = Rw(1,i-L+1);
    Qw(2,i) = Rw(2,i-L+1);
    Qw(3,i) = Rw(3,i-L+1);
end

for i=1:length(Qw)
Q(1,i) = Qw(1,i)/ Qw(3,i);
Q(2,i) = Qw(2,i)/ Qw(3,i);
weights(i) =  Qw(3,i);
end




end

