function v = matvec(M)
% Takes in a matrix and returns a column vector obtained by
% concatenating the columns of the matrix.

s = size(M);
v = zeros(s(1)*s(2),1);

for j = 1:s(2)
    for i = 1:s(1)
        v(s(1)*(j-1) + i) = M(i,j);
    end
end