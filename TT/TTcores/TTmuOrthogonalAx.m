function [Adx] = TTmuOrthogonalAx(A,dU,yl,yr,m,d,r)
%solve A*dx for every search direction,yl and yr are the product before 
%and after each index mu for the TT multiplication Adx
[n_samples,~] = size(A{1});

Adx = zeros(n_samples,d);
for i = 1:d
    dUi = reshape(dU{i},[r(i), m(i), r(i+1)]);
    dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
    AdU = A{i}*dUi;
    AdU = reshape(AdU,n_samples,r(i),r(i+1));

    Adxi = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
    end
    Adx(:,i) = sum(Adxi.*yr{i},2);
end


end

