% ============================================================ %
% For a matrix vec(A) returns A.
% need to specify the number of columns (n2) and rows (n1) in A
% ============================================================ %

function v = unvec(A,n1,n2);
v = [];
for i = 1:n2;
    v(:,i) = A( (i-1)*n1+1 : i*n1,1);
end;
