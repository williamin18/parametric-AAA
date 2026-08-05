function [dXk] = TTcore_Newton(yl,yr,Ak,Xk,residual,lambda)
%TT_NEWTON Summary of this function goes here
%   Detailed explanation goes here



[n_samples,m] = size(Ak);
[~,r_k] = size(yl);
[~,r_k1] = size(yr);

Y_yl = kron( kron(ones(1,r_k1), ones(1,m)), yl );
Y_Ai = kron( kron(ones(1,r_k1), Ak ),       ones(1,r_k));
Y_yr = kron( kron(yr,           ones(1,m)), ones(1,r_k));

Y = Y_yl.*Y_Ai.*Y_yr;

n = r_k*m*r_k1;
%regularization
if lambda ~= 0
    Y = [Y;lambda*eye(n)];
    residual = [residual; -lambda*reshape(Xk,[],1)];
    % residual = [residual; zeros(n,1)];
end


dXk =  Y\residual;
dXk = reshape(dXk,[r_k*m, r_k1]);
end

