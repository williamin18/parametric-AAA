function [g,x2] = TT_Riemannian_fromGTensor(U,V,dUx)
%TT_FORM_RGRAD_TENSOR Summary of this function goes here
%   Detailed explanation goes here
[d,m,r] = TTsizes(U);
g = cell(d,1);

g{1} = [dUx{1} U{1}];
for i = 2:d-1
    temp = zeros(r(i)*2,m(i),r(i+1)*2);
    temp(1:r(i),:,1:r(i+1)) = reshape(V{i}, r(i),m(i),r(i+1));
    temp(r(i)+1:end,:,1:r(i+1)) = reshape(dUx{i}, r(i),m(i),r(i+1));
    temp(r(i)+1:end,:,r(i+1)+1:end) = reshape(U{i}, r(i),m(i),r(i+1));
    g{i} = reshape(temp,r(i)*2*m(i),r(i+1)*2);
end
x2 = g;
g{d} = h2v([v2h(V{d},m(d)); v2h(dUx{d},m(d))],m(d));
x2{d} = h2v([v2h(V{d},m(d)); v2h(dUx{d},m(d))+v2h(U{d},m(d))],m(d));
end

