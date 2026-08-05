function X = i2v(Y,r1)
% i2v  reshape core from inner to vertical unfolding. For a r1*n1*r2
% TT-core, v2z unfold it from  n1*(r1*r2) to (r1*n1)*r2
%
% See also i2v.
[n,rr] = size(Y);
X = reshape(permute(reshape(Y,[n,r1,rr/r1]), [2 1 3]),[r1*n rr/r1]);
end

