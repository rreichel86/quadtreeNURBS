function d = distanceLineSegment(P,A,B)

% unit vector along direction AB
u = (B-A)/norm(B-A);


alpha = (P-A)' * u;

d = norm( (P-A) - alpha*u );


if alpha < 0 
    d = norm(A-P);
end 

if alpha > norm(B-A)
    d = norm(B-P);
end


end 