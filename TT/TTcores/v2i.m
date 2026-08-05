function X = v2i(Y,r1)
% v2i  Reshape core from vertical to inner unfolding. For a r1*n1*r2
% TT-core, v2z unfold it from (r1*n1)*r2 to n1*(r1*r2)
%
% See also i2v.
[nr,r2] = size(Y);
X = reshape(permute(reshape(Y,[r1,nr/r1, r2]), [2 1 3]),[nr/r1 r1*r2]);
end

