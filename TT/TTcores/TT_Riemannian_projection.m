function [Ux] = TT_Riemannian_projection(U,V,x)
%Projection of TT tensors  {U2_1 U2_2 ... U2_k-1 X V2_k+1 ...
%V2_d-1 V2_d} to TT tensors {U_1 U_2 ... U_k-1 Q_lk*X*Q_rk V_k+1 ... V_d-1
%V_d} for every k

[d,m,~] = TTsizes(U);

Y_l = cell(d,1);
Y_l{1} = 1;
Y_r = cell(d,1);
Y_r{d} = 1;

for i = 1:d-1
    gi = h2v( Y_l{i}* v2h(x{i},m(i)) , m(i) );
    Y_l{i+1} = U{i}'*gi;
end

for i = d:-1:2
    gi = x{i}*Y_r{i};
    Y_r{i-1} = v2h(gi,m(i)) * v2h(V{i},m(i))';
end


Ux = cell(d,1);
for i = 1:d-1
    Ux{i} =  h2v( Y_l{i}* v2h(x{i},m(i)) , m(i))* (Y_r{i});
    %Ux{i} =  Ux{i} - U{i}*U{i}'*Ux{i};
end

Ux{d} =  h2v( Y_l{d}* v2h(x{d},m(d)) , m(d));

end

