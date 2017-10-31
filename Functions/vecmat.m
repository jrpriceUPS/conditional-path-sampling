function M = vecmat(v, m, n)
% Takes in a vector and returns an m x n matrix obtained from the
% entries of the given vector.

M = zeros(m,n);

for j = 1:n
    for i = 1:m
        M(i,j) = v(m*(j-1) + i);
    end
end