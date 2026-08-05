function [Ux,V] = TTmuOrthogonalizeRL(U)
%TTMUORTHOGONALIZERL Summary of this function goes here
%   Detailed explanation goes here
[d,m,~] = TTsizes(U);
V = U;
Ux = cell(d,1);
Ux{d} = U{d};
for i = d:-1:2
    [Q, R] = qr(v2h(V{i}, m(i))', 'econ');
    V{i} = h2v(Q', m(i));
    V{i-1} = U{i-1} * R';
    Ux{i-1} = V{i-1};
end

end

