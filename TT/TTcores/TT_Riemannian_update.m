function [x] = TT_Riemannian_update(U,V,dUx,a,max_r)
%TT_RIEMANNIAN_GD Summary of this function goes here
%   Detailed explanation goes here
[d,m,r] = TTsizes(U);
x = cell(d,1);

x{1} = [a*dUx{1} U{1}];
for i = 2:d-1
    temp = zeros(r(i)*2,m(i),r(i+1)*2);
    temp(1:r(i),:,1:r(i+1)) = reshape(V{i}, r(i),m(i),r(i+1));
    temp(r(i)+1:end,:,1:r(i+1)) = reshape(a*dUx{i}, r(i),m(i),r(i+1));
    temp(r(i)+1:end,:,r(i+1)+1:end) = reshape(U{i}, r(i),m(i),r(i+1));
    x{i} = reshape(temp,r(i)*2*m(i),r(i+1)*2);
end
x{d} = h2v([v2h(V{d},m(d)); a*v2h(dUx{d},m(d)) + v2h(U{d},m(d))],m(d));
% x = TTrounding_Randomize_then_Orthogonalize(x,[1 max_r*ones(1,d-1) 1]);
x = TTrounding(x,1e-5,max_r);

end

