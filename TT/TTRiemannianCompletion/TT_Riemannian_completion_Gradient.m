function [V,dUx,dx_TT] = TT_Riemannian_completion_Gradient(A,U,residual,beta,dx_old)
%A is selector, U is current TT-approximation
[n_samples,~] = size(A{1});
[d,m,r] = TTsizes(U);
V = U;
Ux = cell(d,1);
Ux{d} = U{d};


%Right Orthogonalize, find every Ux to compute the partial derivative
V{d} = U{d};
for i = d:-1:2
    [Q, R] = qr(v2h(V{i}, m(i))', 'econ');
    V{i} = h2v(Q', m(i));
    V{i-1} = U{i-1} * R';
    Ux{i-1} = V{i-1};
end

[~,yl] = Ax_left(A,U,d);
[~,yr] = Ax_right(A,V,1);

%Compute Riemannian gradient for each core
dUx = cell(d,1);
for i = 1:d
    for j = 1:m(i)
        temp = residual.*A{i}(:,j);
        temp = yr{i}.*repelem(temp,1,r(i+1));
        dUx{i}((j-1)*r(i)+1:j*r(i),:) = (yl{i}'*temp);
    end
     
end

dx_TT = TT_Riemannian_fromGTensor(U,V,dUx);

%Compute Riemannian gradient vectors evaluated at selected samples
Adx = zeros(n_samples,d);
for i = 1:d
    dUi = reshape(dUx{i},[r(i), m(i), r(i+1)]);
    dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
    AdU = A{i}*dUi;
    AdU = reshape(AdU,n_samples,r(i),r(i+1));

    Adxi = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
    end
    Adx(:,i) = sum(Adxi.*yr{i},2);
end

%Compute step sizes
if beta <=0
    alpha = Adx\residual;
    for i = 1:d
        dUx{i} = alpha(i)*dUx{i};
    end
else
    dU_old = TT_Riemannian_projection(U,V,dx_old);
    Adx_old = zeros(n_samples,d);

    for i = 1:d
        dUi = reshape(dU_old{i},[r(i), m(i), r(i+1)]);
        dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
        AdU = A{i}*dUi;
        AdU = reshape(AdU,n_samples,r(i),r(i+1));

        Adxi = zeros(n_samples,r(i+1));
        for j = 1:r(i+1)
            Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
        end
        Adx_old(:,i) = sum(Adxi.*yr{i},2);
    end
    Adx = [Adx Adx_old];
    alpha = Adx\residual;
    for i = 1:d
        dUx{i} = alpha(i)*dUx{i}+alpha(i+d)*dU_old{i};
    end

end
end

